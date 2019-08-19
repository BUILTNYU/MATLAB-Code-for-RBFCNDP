function [Yset,Yfix,lambda,coeff,Xfix,t,TC]=RBFLOC1final(t0fix,Kfix,alphafix,betafix,A,q,D,B,Ymax,yprime,dquad,Nmax)

%CNDP solution via CORS-RBF developed by Regis & Shoemaker (2007)
%Nmax is the maximum # of function evaluations
%UEAssign is at 100 max iter, .01 percision, m=1 interpolator
%Uses Latin Hypercube Sampling code obtained here: http://www.mathworks.com/matlabcentral/fileexchange/4352 - lhsu.m
%to generate initial set of n0 candidate points
%RBFLOC1 uses multi-start local MSRBF for network without centroids
%Yset and TC are in terms of n, whereas cycle is in the inner loop

start=clock;

weight=0.95;  %only 1 weight
J=size(D,1);
fmax=max(5,J);
rmax=5;

n0=J+1;
n=0;

%create the initial set of n0 candidate points using simple random sampling, with first value defined the same as earlier
TC=0;
%LHS method - i'm aware that it gets skewed a bit
while n<Nmax
    cycle=n0;
    Cfail=0;
    rcount=0;
    ro_n=0.1*Ymax;
    bestTC=inf;
    min_cycle=1;
    n=n+n0;
    if n==n0
        Yset=zeros(J,n0);
    else
        Yset=[Yset zeros(J,n0)];
    end
    s=lhsu(zeros(J,1),Ymax,n0); %this is the LHS code, source cited above
    if n==n0
        Yset=s';
    else
        Yset(:,n-n0+1:n)=s';
    end
    for i=n-n0+1:n
        sumconstraint=(Yset(:,i).^(dquad+1))'*D;
        if sumconstraint>B
            Yset(:,i)=((B/sumconstraint)^(1/(dquad+1)))*Yset(:,i);
        end
    end

    Yfix=zeros(size(A,2),1);

    %Do costly function evaluation
    for i=n-n0+1:n
        for j=1:J
            Yfix(yprime(j,1),1)=Yset(j,i);
        end
        [Xfix,t,TC(i)]=UEAssign2(t0fix,Kfix,alphafix,betafix,Yfix,A,q,100,.01); %in this instance I'm calling a custom UE traffic assignment code to resolve the lower level problem -- obviously this can be swapped out for CNDPs where different lower level problems exist
        TC=[TC;0];
        if TC(i)<bestTC
            min_cycle=i-(n-n0);
            bestTC=TC(i);
        end
    end
     
    while and(rcount<=rmax,n<Nmax)

        %2.1. need to fit/update RBF
        %%let's try using thin-plate splines  
        phi=zeros(cycle,cycle);
        for i=1:cycle
            for j=1:cycle
                if Yset(:,n-cycle+i)-Yset(:,n-cycle+j)==0
                    phi(i,j)=0;
                else
                    phi(i,j)=sum((Yset(:,n-cycle+i)-Yset(:,n-cycle+j)).^2)*log(sqrt(sum((Yset(:,n-cycle+i)-Yset(:,n-cycle+j)).^2))+realmin);  %phi(r)=r^2*log(r)
                end
            end
        end

        %poly is the linear (m=1) polynomial matrix, where poly=[x' 1]
        poly=zeros(cycle,J+1);
        for i=1:cycle
            for j=1:J
                poly(i,j)=Yset(j,n-cycle+i);
            end
            poly(i,J+1)=1;
        end
        lambda=[phi poly; poly' zeros(J+1,J+1)]\[TC(n-cycle+1:n); zeros(J+1,1)];
        coeff=lambda(cycle+1:cycle+J+1);
        lambda=lambda(1:cycle);

        %2.2. then randomly generate candidate points using RBF
        %%use 1000 pts
        candpt=zeros(J,100*J);
        for i=1:100*J
            samples=zeros(J,1);
            for j=1:J
                samples(j,1)=rand;
            end
            candpt(:,i)=Yset(:,n-cycle+min_cycle)+norminv(samples,zeros(J,1),ro_n);
            candpt(:,i)=max(zeros(J,1),candpt(:,i)); %in case it goes less than 0
            sumconstraint=(candpt(:,i).^(dquad+1))'*D;
            if sumconstraint>B
                candpt(:,i)=((B/sumconstraint)^(1/(dquad+1)))*candpt(:,i);
            end
        end

        %2.3. finally, select new function evaluation point
        minval=inf;
        maxval=0;
        mindist=inf;
        maxdist=0;
        dist=zeros(100*J)+inf;
        interpolate=zeros(100*J);
        for i=1:100*J
            for j=1:cycle
                interpolate(i)=interpolate(i)+lambda(j)*sum((candpt(:,i)-Yset(:,n-cycle+j)).^2)*log(sqrt(sum((candpt(:,i)-Yset(:,n-cycle+j)).^2))+realmin);
                %determining distance
                dist(i)=min(dist(i),sqrt(sum((candpt(:,i)-Yset(:,n-cycle+j)).^2)));            
            end
            interpolate(i)=interpolate(i)+candpt(:,i)'*coeff(1:J)+coeff(J+1);
            if interpolate(i)<minval
                minval=interpolate(i);     %set the minimum value
            end

            if interpolate(i)>maxval
                maxval=interpolate(i);     %set the max value and i
            end

            if dist(i)<mindist
                mindist=dist(i);
            end

            if dist(i)>maxdist
                maxdist=dist(i);
            end
        end
        score=inf;
        min_i=1;
        for i=1:1000
            if maxval~=minval
                Vr=(interpolate(i)-minval)/(maxval-minval);
            else
                Vr=1;
            end
            if maxdist~=mindist
                Vd=(maxdist-dist(i))/(maxdist-mindist);
            else
                Vd=1;
            end

            if weight*Vr+(1-weight)*Vd<score
                score=weight*Vr+(1-weight)*Vd;
                min_i=i
            end
        end

        %2.4 & 2.5 Evaluate costly function and update information

        Yset=[Yset candpt(:,min_i)];
        Yfix=zeros(size(A,2),1);
        for i=1:J
            Yfix(yprime(i,1),1)=Yset(i,n+1);
        end
        [Xfix,t,TC(n+1)]=UEAssign2(t0fix,Kfix,alphafix,betafix,Yfix,A,q,100,.01); %here too - custom UE traffic assignment code is called
        if TC(n+1)<bestTC
            bestTC=TC(n+1)
            min_cycle=cycle+1;
            Cfail=0;
        else
            Cfail=Cfail+1;
        end
        if Cfail>fmax
            Cfail=0;
            rcount=rcount+1;
            ro_n=0.5*ro_n;
        end
        if n+1<Nmax
            TC=[TC;0];
        end
        n=n+1
        cycle=cycle+1
    end
end


%3. the best solution is obtained at Nmax.


function [f,Cn,fc]=densityWatsonKernel(XYZ,w,grid,r,M0,ifig);
% ifig=4;

%% Find the self-consistent smoothing parameter
nPt=size(XYZ,2);
if nPt>1000;
    % %%
    % sIterKeep=[10,20,50,100,200,500,1000];
    % CnKeep=NaN(size(sIterKeep));
    % for iN=1:numel(sIterKeep);
    %     sIter=sIterKeep(iN)
    % repeat nIter times with subset of size sIter
    nIter=10; CnIter=NaN(1,nIter); sIter=100;
    for iIter=1:nIter;
        iSelect=ceil(nPt*rand(1,sIter));
        CnIter(iIter)=smoothingSelfConsistent(XYZ(:,iSelect),w(iSelect),ifig);
    end
    % CnKeep(iN)=mean(CnIter);
    Cn=round(mean(CnIter)*(nPt/sIter)^(1/3),1,'significant');
    % end;
    %
else
    Cn=round(smoothingSelfConsistent(XYZ,w,ifig),1,'significant');
end

%%


% Cn=smoothingSelfConsistent(XYZ,ifig);
% estimate the density; no weighing
f=densityEstimate(XYZ,w,grid.XYZ,Cn);
% Find the contour level of the density estimated
C=decideContour(grid.XYZ,f,r,M0,ifig);%
% identify where the contours are in the projected Schmidt plane
% [M,H]=contourf(grid.l,grid.l,reshape(f,size(grid.x)),C);
M=contourc(grid.l,grid.l,reshape(f,size(grid.x)),C);
Mn=M;
fc=[];ic=0;nM=size(Mn,2);
while nM>=2; % there is at least one more contour;
    ic=ic+1;nn=Mn(2,1);
    fc(ic).level=Mn(1,1);
    fc(ic).number=nn;
    fc(ic).contour.x=Mn(1,1+[1:nn]);
    fc(ic).contour.y=Mn(2,1+[1:nn]);
    fc(ic).contour=schmidt2sphere(fc(ic).contour);
    fc(ic).contour.lower=truncateLower(fc(ic).contour);
    Mn=Mn(:,[nn+2:end]);
    nM=size(Mn,2);
end
end
function Cn=smoothingSelfConsistent(P3,w,ifig);
n=size(P3,2);
cosTheta=(P3'*P3); %xi*xj+yi*yj+zi*zj

figure(ifig); subplot 121; hold on;
Cr=10;
Cmin=1;Cmax=Cr*Cmin;nC=10;CV=false;
while ~CV
    CAll=10.^(linspace(log10(Cmin),log10(Cmax),nC));
    Lstore=NaN(size(CAll));
    for iC=1:numel(CAll);
        Lstore(iC)=pseudoLogLikelihood(cosTheta,w,CAll(iC),n);
    end
    plot(CAll,Lstore,'.-')
    % Update the estimates of C
    [Lmax,iLmax]=max(Lstore);
    if iLmax==1;
        Cmax=Cmin*Cr/2;
        Cmin=Cmin/Cr;
    elseif iLmax==nC;
        Cmin=Cmax/(Cr/2);
        Cmax=Cmax*Cr;
    else
        if 2*Lmax/(Lstore(iLmax+1)+Lstore(iLmax-1))>=0.99;
            CV=true;
        else
            Cmin=CAll(iLmax-1);
            Cmax=CAll(iLmax+1);
        end
    end
    if or(Cmin<1e-2,Cmax>1e3);
        CV=true;
    end
end
Cn=(CAll(iLmax));

plot(Cn,Lmax,'ok');
box on; xlabel('Smoothing parameter C_n'); ylabel('pseudo-log-likelihood L(C_n)')


    function L=pseudoLogLikelihood(P,w,C,n);
        % W2=(C/(4*pi*n*sinh(C)))*exp(C*P);
        if C>0;
            Cw=2*pi^(3/2)*erfi(sqrt(C))/sqrt(C);
        else
            Cw=2*pi^(3/2)*erf(sqrt(-C))/sqrt(-C);
        end
        % W2=exp(C*(P.^2))/(n*Cw);
        W2=repmat(w,size(P,1),1).*exp(C*(P.^2))/Cw;
        fIncomplete=NaN(1,n);
        for i=1:n;
            iInc=[1:i-1,i+1:n];
            fIncomplete(i)=sum(W2(i,iInc))/sum(w(iInc));
        end
        L=sum(log(fIncomplete));
    end
end

%%
function f=densityEstimate(XYZ,w,XYZi,C);
n=size(XYZ,2); %size of the dataset
ng=size(XYZi,2); %size of the grid
wt=sum(w); %total weight

D2=(XYZi'*XYZ); %xi*xj+yi*yj+zi*zj
% W2=(C/(4*pi*n*sinh(C)))*exp(C*D2);
if C>0;
    Cw=2*pi^(3/2)*erfi(sqrt(C))/sqrt(C);
else
    Cw=2*pi^(3/2)*erf(sqrt(-C))/sqrt(-C);
end
W2=repmat(w,ng,1).*exp(C*(D2.^2))/Cw;
f=sum(W2,2)/wt;
end

%%
function C=decideContour(XYZ,f,r,M0,ifig);
ik=0;
M=M0;
iLower=find(XYZ(3,:)<0);
fsort=sort(f(find(XYZ(3,:)<0)));
N=numel(fsort);
jn=N;j=[];ij=0;
while jn>M;
    ij=ij+1;
    jn=jn-M;
    j(ij)=jn;
    M=M*r;
end
j=round(j);

figure(ifig); subplot 122; hold on;
plot([1:N]/N,fsort);
plot(j/N,fsort(j),'o');
box on; xlabel('Fraction of data'),ylabel('density estimate');
title(sprintf('r=%g, M_0=%g',r,M0));

C=fsort(j);
end

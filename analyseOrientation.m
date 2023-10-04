function orientation=analyseOrientation(X,Y,Z,W)
% orientation=analyseOrientation(X,Y,Z,W)
% input: data set
%   X, Y, Z: coordinates on unit sphere
%   W: weight for each data point
% output: orientation structure
%   orientation.Matrix
%   orientation.eigenValue
%   orientation.eigenVectors
%   orientation.eigenAzimuth in radians
%   orientation.eigenElevation in radians, 0 at equator
%   orientation.eigenRadius should be ones
%   orientation.shapeFactor=log(E(3)/E(2))/log(E(2)/E(1));
%   orientation.shapeStrength=log(E(3)/E(1));
%   orientation.cone3: 3-D uncertainty cone around highest eigenvector.
%   orientation.cone1: 3-D uncertainty cone around smallest eigenvector.


% Construct the orientation Matrix (weighted average)
orientationMatrix=zeros(3);
orientationMatrix(1,1)=sum(X.*X.*W);
orientationMatrix(1,2)=sum(X.*Y.*W);
orientationMatrix(1,3)=sum(X.*Z.*W);
orientationMatrix(2,1)=sum(Y.*X.*W);
orientationMatrix(2,2)=sum(Y.*Y.*W);
orientationMatrix(2,3)=sum(Y.*Z.*W);
orientationMatrix(3,1)=sum(Z.*X.*W);
orientationMatrix(3,2)=sum(Z.*Y.*W);
orientationMatrix(3,3)=sum(Z.*Z.*W);
Wtotal=sum(W);
orientationMatrix=orientationMatrix/Wtotal;
% eigenvalue and eigenvectors
[eigVector,eigDiagonal]=eig(orientationMatrix);
eigValue=diag(eigDiagonal);
[E,iE]=sort(eigValue);
eigVector=eigVector(:,iE);
eigValue=E;
eigDiagonal=diag(E);

% Flip eigenvectors so they always point upward/
for idown = find(eigVector(3,:)<0)
    eigVector(:,idown)=-(eigVector(:,idown));
end
% reproject 
for i=1:3
    [eigAzimuth(i),eigElevation(i),eigRadius(i)]=cart2sph(eigVector(1,i),eigVector(2,i),eigVector(3,i));
end
% store
orientation.Matrix=orientationMatrix;
orientation.eigenValue=eigValue;
orientation.eigenVectors=eigVector;
orientation.eigenAzimuth=eigAzimuth;
orientation.eigenElevation=eigElevation;
orientation.eigenRadius=eigRadius;

% E=sort(eigValue);
orientation.shapeFactor=log(E(3)/E(2))/log(E(2)/E(1));
orientation.shapeStrength=log(E(3)/E(1));

% R = sqrt((sum(X.*W)).^2+(sum(Y.*W)).^2+(sum(Z.*W)).^2)./sum(W);+

[orientation.cone3,orientation.cone1]=interpretOrientation(numel(W),eigValue,eigVector,[X;Y;Z],W);

end

function [cone3,cone1]=interpretOrientation(n,eigValue,eigVectors,XYZ,W)

% Any axis vs. random
% n=size(XYZ,2);
Sig=[10,5,2.5,1]; %significance in table 
Bplus=[0.788,0.874,0.948,1.038]; % b coefficient for each significance
b=(eigValue-(1/3))*sqrt(n); % scale with DOF
%% Bipolar distribution Anderson and Stephens (1972)
bmax=b(3);
if bmax>max(Bplus)
    %disp(sprintf('Random over bipolar (any axis): less than %g %%',min(Sig)));
elseif bmax<min(Bplus)
    %disp(sprintf('Random over bipolar (any axis): more than %g %%',max(Sig)));
else
    %disp(sprintf('Random over bipolar (any axis): %g %%',interp1(Bplus,Sig,bmax)));
end
%% Girdle distribution
bmin=-b(1);
if bmin>max(Bplus)
    %disp(sprintf('Random over girdle (any axis): less than %g %%',min(Sig)));
elseif bmin<min(Bplus)
    %disp(sprintf('Random over girdle (any axis): more than %g %%',max(Sig)));
else
    %disp(sprintf('Random over girdle (any axis): %g %%',interp1(Bplus,Sig,bmin)));
end
%%
V1=eigVectors(:,1);
V2=eigVectors(:,2);
V3=eigVectors(:,3);
Wt=sum(W);

%% random vs. specific eigenvectors
% Sn3=sum((V3'*XYZ).^2)/n; %Scale with weight
Sn3=sum(W.*(V3'*XYZ).^2)/Wt; %Scale with weight
S3=sqrt(n)*(Sn3-1/3)/sqrt(4/45); %scale with DOF
%disp(sprintf('Random over bipolar at the 3rd eigenvector: %g%% (S3 = %g)',100*(1/2)*erfc(S3/sqrt(2)),S3));

% Sn1=sum((V1'*XYZ).^2)/n; %Scale with weight
Sn1=sum(W.*(V1'*XYZ).^2)/Wt; %Scale with weight
S1=sqrt(n)*(Sn1-1/3)/sqrt(4/45); %scale with DOF
%disp(sprintf('Random over bipolar at the 1rd eigenvector: %g%% (S1 = %g)',100*(1-(1/2)*erfc(S1/sqrt(2))),S1));

%
% Sn=sum((V3'*XYZ).^2)/n;S=n^(1/3)*(Sn-1/3)/(4/45)^(1/2);X=erfc(S/sqrt(2))/2;
% Sn=sum(W.*(V3'*XYZ).^2)/Wt;S=n^(1/2)*(Sn-1/3)/(4/45)^(1/2);X=erfc(S/sqrt(2))/2;
% %disp(sprintf('Probability of random over 3rd eigenvector %g %%',round(X*100,2,'significant')));
%%
e11=sum(W.*((V1'*XYZ).^2).*((V3'*XYZ).^2))/((eigValue(1)-eigValue(3))^2);
e22=sum(W.*((V2'*XYZ).^2).*((V3'*XYZ).^2))/((eigValue(2)-eigValue(3))^2);
e12=sum(W.*(V1'*XYZ).*(V2'*XYZ).*((V3'*XYZ).*2))/((eigValue(1)-eigValue(3))*(eigValue(2)-eigValue(3)));
E=[e11,e12;e12,e22]/Wt;
% F=inv(E);
alp=0.1;
cone3=makeCone(E,n,alp,eigVectors);

e11=sum(W.*((V1'*XYZ).^2).*((V3'*XYZ).^2))/((eigValue(1)-eigValue(3))^2);
e22=sum(W.*((V2'*XYZ).^2).*((V1'*XYZ).^2))/((eigValue(2)-eigValue(1))^2);
e12=sum(W.*(V3'*XYZ).*(V2'*XYZ).*((V1'*XYZ).*2))/((eigValue(1)-eigValue(2))*(eigValue(1)-eigValue(3)));
E=[e11,e12;e12,e22]/Wt;
% F=inv(E);
alp=0.1;
cone1=makeCone(E,n,alp,fliplr(eigVectors));

%% Significance of principal axes
Fstar=zeros(3); Fstar([2:3],[2:3])=inv(cone3.E);
H=eigVectors;
X2=V3'*H*Fstar*transpose(H)*V3*n;
%disp(sprintf('Significance probability that 3rd eigenvector is the principal axis: %g%%',100*exp(-X2/2)));

Fstar=zeros(3); Fstar([1:2],[1:2])=inv(cone1.E);
H=fliplr(eigVectors);
X2=V1'*H*Fstar*transpose(H)*V1*n;
%disp(sprintf('Significance probability that 1rd eigenvector is the polar axis: %g%%',100*exp(-X2/2)));

%% Parameter estimation
Df=@(z)exp(z)./(sqrt(pi*z).*erfi(sqrt(z)))-1./(2*z); %6.21
KappaBipolar=fzero(@(x)Df(x)-eigValue(3),1);
%disp(sprintf('Bipolar parameter estimate: %g',KappaBipolar));

KappaGirdle=fzero(@(x)Df(x)-eigValue(1),-1);
%disp(sprintf('Girdle parameter estimate: %g',KappaGirdle));

cosSq1=(V1'*XYZ).^2;
cosSq3=(V3'*XYZ).^2;
figure(4); clf; 
subplot 121;
plot(log(n./(n+(1/2)-[1:n])),sort(1-cosSq3)); hold on; plot([0,KappaBipolar],[0,1]); 
title(sprintf('Q-Q plot, bipolar axis, Kappa=%g',KappaBipolar));
xlabel('Test quartile'); ylabel('Sample quartile');
subplot 122;
plot(chi2inv(([1:n]-(1/2))/n,1),sort(cosSq1)); hold on; plot([0,abs(2*KappaGirdle)],[0,1]); 
title(sprintf('Q-Q plot, girdle pole, Kappa=%g',KappaGirdle));
xlabel('\chi_1^2 quartile'); ylabel('Sample quartile');
%%
% Fstar=[[0,0,0];[[0;0],F]];
% X2=n*V3'*eigVectors*Fstar*eigVectors'*V3;
% A2=exp(-X2/2);
% %disp(sprintf('significance probability that V3 is the true axis: %g',round(A2,3,'significant')));
%%
% random vs. bipolar
% bmin=min(b);
% Bcrit=-Bcrit;
% if bmin<min(Bcrit);
%     %disp('Random over girdle: less than 1%')
% elseif bmin>max(Bcrit);
%     %disp('Random over girdle: more than 10%')
% else
%     fprintf('Random over girdle: %g%%\n',round(interp1(Bcrit,P,bmin),2,'significant'));
% end

%%
    function cone=makeCone(E,n,alp,eigVectors);
        % Build Z matrix 3.21
        F=inv(E);
        A=F(1,1); B=F(1,2); C=F(2,2); D=-2*log(alp)/n;
        Z=[A,B;B,C];
        % eigenvectors 3.22, 3.23
        [Y,T]=eig(Z);
        % sin of vertical angles
        g1=sqrt(D/T(1,1));
        g2=sqrt(D/T(2,2));
        % angles for a nice graph
        psi=linspace(0,2*pi,60);
        % ellipse centered on the pole
        v1=g1*cos(psi); v2=g2*sin(psi);
        x=Y(1,1)*v1-Y(2,1)*v2;
        y=Y(2,1)*v1+Y(1,1)*v2;
        z=sqrt(1-x.^2-y.^2);
        if isreal(z);
            % rotate to center on eigenvectors
            xyzr=eigVectors*[x;y;z];
            % move to lower hemisphere;
            % break the line when crossing the equator.
            iEquator=find(diff(sign(xyzr(3,:)))~=0);
            if ~isempty(iEquator)
                for iCross=fliplr(iEquator);
                    fr=-xyzr(3,iCross)/(xyzr(3,iCross+1)-xyzr(3,iCross));
                    xc=xyzr(1,iCross)+fr*(xyzr(1,iCross+1)-xyzr(1,iCross));
                    yc=xyzr(2,iCross)+fr*(xyzr(2,iCross+1)-xyzr(2,iCross));
                    xyzr=[xyzr(:,1:iCross),...
                        sign(xyzr(3,iCross+1)-xyzr(3,iCross))*[[xc,NaN,-xc];[yc,NaN,-yc];[0,NaN,0]],...
                        xyzr(:,iCross+1:end)];
                end
            end
            % move upper hemisphere points to the lower hemisphere.
            iUpper=find(xyzr(3,:)>0);
            xyzr(:,iUpper)=-xyzr(:,iUpper);
            % Project to Schmidt
            colat=acosd(xyzr(3,:));
            lon=atan2d(xyzr(2,:),xyzr(1,:));
            [xs,ys]=sphere2schmidt(colat,lon);
            %disp('What do you mean, z is not real?')
        else % ellipse is too big; use the equator
            xc=cos(psi);yc=sin(psi);zc=sin(psi)*0;
            xyzr=[xc,yc,zc];
            lon=psi;colat=psi*0+pi/2;
            xs=xc*sqrt(2);ys=yc*sqrt(2);
        end

        % record for output
        cone.XYZ=xyzr;
        cone.lon=lon;cone.colat=colat;
        cone.x=xs;cone.y=ys;
        cone.verticalAngles=[asind(g1),asind(g2)];
        cone.E=E;

    end
end

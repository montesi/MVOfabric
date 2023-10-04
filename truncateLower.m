function trunc=truncateLower(cont);
colat=cont.colat;
iLow=find(colat>90);
colatT=colat(iLow);
lonT=cont.lon(iLow);
rho = 2*sind((180-colatT)/2);


iUpper=find(colat<90);
colatT=colat; colatT(iUpper)=90;
lonT=cont.lon;
rho = 2*sind((180-colatT)/2);

%Close contour along sphere's edge
if and(colatT(end)==90,colatT(1)==90)
    dLon=lonT(1)-lonT(end);
    nloopBack=ceil(dLon/5);
    loopBackLon=linspace(lonT(end),lonT(1));
    loopBackColat=90+loopBackLon*0;
    loopBackRho=sqrt(2)+loopBackLon*0;
    lonT=[lonT,loopBackLon];
    colatT=[colatT,loopBackColat];
    rho=[rho,loopBackRho];
end


trunc.x = rho.*cosd(lonT);
trunc.y = rho.*sind(lonT);
trunc.lon=lonT;
trunc.colat=colatT;


trunc.XYZ=[sind(colatT).*cosd(lonT);...
    sind(colatT).*sind(lonT);...
    cosd(colatT)];
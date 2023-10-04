function [x,y] = sphere2schmidt(colat,lon);

rho = 2*sind((180-colat)/2);
x = rho.*cosd(lon);
y = rho.*sind(lon);
function grid=schmidt2sphere(grid);
% convert a set of points in a Schmidt project to their 3D projection on a
% unit sphere
x=grid.x;y=grid.y;

rho=(x.^2+y.^2).^(1/2);

lon=atan2d(y,x);
colat=180-2*asind(rho/2);
grid.lon=lon;grid.colat=colat;

X=sind(colat).*cosd(lon);
Y=sind(colat).*sind(lon);
Z=cosd(colat);

grid.XYZ=[X(:),Y(:),Z(:)]';

function grid=SchmidtProjection(data,ifig,iDet,dAngle,nx,OR);%,dataPerc,BDW)
% SchmidtProjection(data,ifig,dAngle,nx);
% 1) Project axial data on a Schmidt net;
% 2) Define contour data with self-consistent concentration parameter and
% both Watson and Fisher kernels. 
% Current setup is Watson with weights, Watson without weights
% Inputs: 
%   data.Phi, data.Theta: colatitude and longitude of axial data
%   ifig: Figure number for graphical output
%   iDet: Figure number for the determination of smoothing parameters and contours
%   dAngle: spacing of reference grid for the Schmidt net
%   nx: number of sample in x (and y) points across the Schmidt net.
% Ouput: grid
%   grid.l: vector in x or y directions of grid point location
%   grid.x, grid.y: grid points on Schmidt net
%   grid.lon, grid.colat: grid points longitude and colatitude
%   grid.XYZ: grid points coordinates on a unit sphere (3 x nx^2 matrix)
%   grid.fWatson: density esitmate using a Watson kernel
%   grid.CnWatson: Concentration factor for the density estimate using a Watson kernel
%   grid.fcWatson: Contours of the Watson density estimate (1xnumber of contours structure)
%       grid.fcWatson.level: level of a particular contour      
%       grid.fcWatson.number: number of points a particular contour
%       grid.fcWatson.contour: information of a particular
%           grid.fcWatson.contour.x, grid.fcWatson.contour.y: points on Schmidt net
%           grid.fcWatson.contour.lon, grid.fcWatson.contour.colat: longitude and colatitude 
%           grid.fcWatson.contour.XYZ: coordinates in the unit sphere
%           grid.fcWatson.contour.lower: information about the contour restricted to the lower hemisphere
%               grid.fcWatson.contour.lower.x, grid.fcWatson.contour.lower.y: points on Schmidt net
%               grid.fcWatson.contour.lower.lon, grid.fcWatson.contour.lower.colat: longitude and colatitude 
%               grid.fcWatson.contour.lower.XYZ: coordinates in the unit sphere
%   grid.*WatsonNoWeight or grid.*Fisher alternative versions of the grid results

% Project points on a Schmidt net
[x,y] = sphere2schmidt(data.Phi,data.Theta);

% Make grid in projected 2D plane first 
li=linspace(-1,1,nx)*sqrt(2);
[xi,yi]=meshgrid(linspace(-1,1,nx)*sqrt(2));
grid.l=li;grid.x=xi;grid.y=yi;
% Find where these points are on the unit sphere
grid=schmidt2sphere(grid);

% initialize figure
figure(ifig); clf; hold on;
colormap sky;

ncase=numel(data.case);
nsub=ncase+2;
nsx=floor(sqrt(nsub));
nsy=ceil(nsub/nsx);


%% Plot the axial data
subplot (nsy,nsx,1); hold on;
plot(x,y,'.');
axis equal;
plotAddSchmidtGrid(gca,dAngle); %reference frame
if isfield(data,'Label');
text(0,1.7,data.Label,...
    'FontWeight','bold',...
    'HorizontalAlignment','center',...
    'fontsize',12)
else
text(0,1.7,'All vectors',...
    'FontWeight','bold',...
    'HorizontalAlignment','center',...
    'fontsize',12)
end
set(gca,'Visible','off')
axis tight;
colorbar

%% Density contours, Fisher kernel
% [f,Cn,fc]=densityFisherKernel(data.XYZ,grid,2,10,4);

%% Density contours, Watson kernel
% [f,Cn,fc]=densityWatsonKernel(data.XYZ,ones(size(data.Phi)),grid,1.5,10,2);
figure(iDet); clf; 
for icase=1:ncase;
    [f,Cn,fc]=densityWatsonKernel(data.XYZ,data.case(icase).weight,grid,1.5,10,iDet);
    grid.case(icase).fWatson=f;grid.case(icase).fcWatson=fc;grid.case(icase).CnWatson=Cn;
    % plot the contours
    figure(ifig); subplot (nsy,nsx,icase+1); hold on;
    SchmidtGeneral(f,Cn,fc,data.case(icase).label,dAngle,data.case(icase).orientation);
axis tight;
end
figure(ifig); subplot (nsy,nsx,nsub); hold on;
plotOrientation(OR');
end

function SchmidtGeneral(f,Cn,fc,label,dAngle,orientation)
% f: density estimate on grid
% Cn: Concentration factor for the density estimate
% fc: Contours of the density estimate (1xnumber of contours structure)
% label: text for the title
% dAngle: controls the graticule
% orientation: direction of the 3 eigenvectors

% plot the contours
colMap=colormap(sky);
fMax=max(f(:)); 
colInd=linspace(0,fMax,size(colMap,1));
for ic=1:numel(fc); 
    colNow=interp1(colInd,colMap,fc(ic).level);
    plot(fc(ic).contour.lower.x,fc(ic).contour.lower.y,'color',colNow)
end
axis equal;
plotAddSchmidtGrid(gca,dAngle);  %reference frame
if ~isempty(orientation);
    addOrientation(orientation,gca);
end
text(0,1.7,sprintf('%s, C=%g',label,Cn),...
    'FontWeight','bold',...
    'HorizontalAlignment','center',...
    'fontsize',12)
set(gca,'Visible','off','clim',[0,fMax]);
% set(gca,'Clim',[0,fmax]);
colorbar('Ticks',flipud(unique([fc.level])));
end

function addOrientation(orientation,hPlot)
% Add eigenvectors and cone to a graph
% input: 
%   orient: orientation structures
%       orient.eigenElevation
%       orient.eigenAzimuth
%       orient.cone3: uncertainty cone, already projected
%           orient.x, orient.y: Schmidt net coordinates%           
%   hPlot: axis handle for plotting
% NEED sphere2schmidt

% find colatitude and longitude
eigColat=90-rad2deg(orientation.eigenElevation);
eigLon=rad2deg(orientation.eigenAzimuth);
% move to lower sphere
iTop=find(eigColat<90);
eigColat(iTop)=180-eigColat(iTop);
eigLon(iTop)=mod(180+eigLon(iTop),360);
% Project on Schmidt net
[X,Y]=sphere2schmidt(eigColat,eigLon);
% Include axes
scatter(X,Y,orientation.eigenValue*500,'or','filled','Parent',hPlot,'MarkerFaceAlpha',0.5);
% Add uncertainty cone
if ~isempty(orientation.cone3)
    plot(orientation.cone3.x,...
        orientation.cone3.y,...
        'r','linewidth',2)
end
if ~isempty(orientation.cone1)
    plot(orientation.cone1.x,...
        orientation.cone1.y,...
        'r','linewidth',2)
end
end

%%
function plotAddSchmidtGrid(ax,dAngle)
% Add reference lines on Schmidt net
% Inputs
%       ax: axis handle
%       dAngle: spacing of reference lines, in degree

% First, spokes
nSpokes=ceil(360/dAngle)+1;
thSpokes=linspace(0,360,nSpokes);
stretch=1.05;
for iSpokes=1:nSpokes-1;
    thNow=thSpokes(iSpokes);
    if cosd(thNow)>0; HA='left';
    elseif cosd(thNow)<0; HA='right';
    else HA='center';
    end
    plot([0,sqrt(2)]*cosd(thNow),...
        [0,sqrt(2)]*sind(thNow),...
        'k','linewidth',0.5)
    text(stretch*sqrt(2)*cosd(thNow),...
        stretch*sqrt(2)*sind(thNow),...
        sprintf('%i%s',thNow,'\circ'),...num2str(thNow),...%'fontsize',8,...
        'HorizontalAlignment',HA)
end

% Second, circles
nCircles=ceil(90/dAngle)+1;
phiCircles=linspace(0,90,nCircles);
rhoCircles=2*cosd((180-phiCircles)/2);
Angles=linspace(0,360,100);
thNow=(thSpokes(1)+thSpokes(2))/2;
for iCircle=1:nCircles
    rhoNow=rhoCircles(iCircle);
    plot(rhoNow*cosd(Angles),...
        rhoNow*sind(Angles),...
        'k','linewidth',0.5)
    phAng=round(phiCircles(iCircle),0); %display colatitude
    if cosd(thNow)>0; HA='left';
    elseif cosd(thNow)<0; HA='right';
    else HA='center';
    end
    text(stretch*rhoNow*cosd(thNow),stretch*rhoNow*sind(thNow),...
        sprintf('%i%s',phAng,'\circ'),...num2str(thNow),...%'fontsize',8,...
        'HorizontalAlignment',HA,'VerticalAlignment','bottom')
end
end

function plotOrientation(OR);
% plotOrientation(OR,ifig)
% show shape factor and intensity on axial data discriment plot
% Inputs: 
%   OR: Structure including the eigenvalues of data to plot
%       OR is a structure, one element for each series to  plot
%       OR.eigVal: 3x1 vector with eigenvalues
%       OR.shape: Symbol to use
%       OR.color: color to use
%       OR.label: series label
%   ifig: Figure number to use. Figure is completely overwritten

pMax=2; % axis limit;
% figure(ifig); 
% clf; hold on;
% shape factors
for gamma=[0,0.2,0.5,1,2,5];
    ang=atan(gamma);
    plot([0,cos(ang)]*pMax,[0,sin(ang)]*pMax,'k');
    text(cos(ang)*pMax,sin(ang)*pMax,...
        sprintf('%s\\gamma = %g','  ',gamma),...
        'Rotation',rad2deg(ang),'VerticalAlignment','baseline');
end
% fabric intensity
for zeta=0:0.5:pMax
    plot([0,zeta],[zeta,0],'k')
    text(zeta/2,zeta/2,...
        sprintf('%s\\zeta = %g','  ',zeta),...
        'Rotation',-45,'VerticalAlignment','bottom');
end
xlabel('log(\tau_2/\tau_1)')
ylabel('log(\tau_3/\tau_2)')
axis equal; box on;
axis([0,pMax,0,pMax]*1.2)
% data to plot
nOR=numel(OR);
Hp=NaN(size(OR));
for iOR=1:nOR
    if ~isempty(OR(iOR).eigVal)
    Hp(iOR)=plot(log(OR(iOR).eigVal(2)/OR(iOR).eigVal(1)),...
        log(OR(iOR).eigVal(3)/OR(iOR).eigVal(2)),...
        'MarkerSize',15,'Color','k','linestyle','none',...
        'marker',OR(iOR).shape,'MarkerFaceColor',OR(iOR).color,...
        'DisplayName',OR(iOR).label);
    end
end
legend(Hp(find(~isnan(Hp))),'Location','northeast');
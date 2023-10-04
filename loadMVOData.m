function [dataSub,dataPore] = loadMVOData(fName);
% (dataSub, dataPore) = loadMVOData('200cubes_POREINERTIA.csv')
% Define the name of each pore-related column
% PORESHAPECOL = 'PoreShape' % name of the column that includes each pore's shape parameter
% POREVOLCOL = 'PoreVolume' % name of the column that includes each pore volume
% POREPHICOL = 'Column_1' % name of the column that includes the latitude for each pore volume
% PORETHETACOL = 'Column_0' % name of the column that includes the longitude for each pore volume

% Load PERGEOS CSV output
A=readtable(fName);
%% Parse and store subvolume information
numSub=size(A,1);
dataSub.meltFraction = A.TotalMelt;
dataSub.xConnected = A.XConnectedMelt;
dataSub.yConnected = A.YConnectedMelt;
dataSub.zConnected = A.ZConnectedMelt;
dataSub.xPermeability = A.kX;
dataSub.yPermeability = A.kY;
dataSub.zPermeability = A.kZ;
dataSub.xPosition = A.Xi;
dataSub.yPosition = A.Yi;
dataSub.zPosition = A.Zi;

%% Parse pore information
if max(strcmp(A.Properties.VariableNames,'PoreVolume')); %in case it's not recorded.
    % Initialize information to keep
    volPhi = []; % Pore orientation colatitude
    volTheta = []; % Pore orientation longitude
    volShape = []; % Pore shape parameter
    volVol = []; % Pore volume
    volSub = []; % Subvolume in which the pore was defined.
    % Compile information
    for iSub = 1:numSub
        nvol = numel(str2num(cell2mat(A.PoreVolume(iSub))));
        volPhi=cat(2,volPhi,str2num(cell2mat(A.Column_1(iSub))));
        volTheta=cat(2,volTheta,str2num(cell2mat(A.Column_0(iSub))));
        volShape=cat(2,volShape,str2num(cell2mat(A.PoreShape(iSub))));
        volVol=cat(2,volVol,str2num(cell2mat(A.PoreVolume(iSub))));
        volSub=cat(2,volSub,iSub*ones(1,nvol));
    end

    % 3D position
    volX=sind(volPhi).*cosd(volTheta);
    volY=sind(volPhi).*sind(volTheta);
    volZ=cosd(volPhi);
    volXYZ=[volX;volY;volZ];
    % move to lower sphere
    iUpper = find(volZ>0);
    volXYZ(:,iUpper)=-volXYZ(:,iUpper);

    dataPore.XYZ=volXYZ;
    dataPore.Phi=acosd(volXYZ(3,:));
    dataPore.Theta=atan2d(volXYZ(2,:),volXYZ(1,:));


    dataPore.Volume = volVol;
    dataPore.Shape = volShape;
    % dataPore.meltPhi = 90-volPhi;
    % dataPore.Phi = volPhi; %colatitude
    % dataPore.Theta = mod(volTheta,360); %longitude
    dataPore.Sub = volSub;
    % dataPore.XYZ = volXYZ;
    % dataPore.X=sind(dataPore.meltPhi).*cosd(volTheta);
    % dataPore.Y=sind(dataPore.meltPhi).*sind(volTheta);
    % dataPore.Z=cosd(dataPore.meltPhi);
else
    dataPore=[];
end

end
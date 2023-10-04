% % Need a distribution to start with
% RepoInfo=dir('../osfstorage-archive_QC0705/*.csv');
% nRepo=numel(RepoInfo);
% for iRepo=1:nRepo
%     [data(iRepo).Sub, data(iRepo).Pore] = loadMVOData([RepoInfo(iRepo).folder,'/',...
%         RepoInfo(iRepo).name]);
%     data(iRepo).Name=RepoInfo(iRepo).name;
% end
% save('CUBES.mat','data');
load CUBES.mat; nRepo=numel(data);

%% orientation
clear OR
for iRepo=1:nRepo;
    if ~isempty(data(iRepo).Pore)
        [F,N,E]=fileparts(data(iRepo).Name);
        disp(sprintf('\n\nWorking on %s',N));
        % Orientation data label
        data(iRepo).Pore.Label=N;
        [data(iRepo).Pore,ORn]=analyseSequence(data(iRepo).Pore,iRepo/nRepo);
        for icase=1:numel(ORn)
            OR(icase,iRepo)=ORn(icase);%OR(2,iRepo)=ORn(2);OR(3,iRepo)=ORn(3);
        end
    end
end
%%
figure(5); clf; hold on;
plotOrientation(OR');
print(5,'-dpdf','Shape_AllData_2','-bestfit');
disp('All done!')
return
%%
function [data,OR]=analyseSequence(data,bright);
% Define different weight cases
% data.case(1).weight=ones(size(data.Phi)); data.case(1).label='Unweighted data';data.case(1).case='Unweighted';
data.case(1).weight=data.Volume; data.case(1).label='Data weighted by volume';data.case(1).case='Volume';
% data.case(3).weight=data.Shape; data.case(3).label='Data weighted by shape';data.case(3).case='Shape';

% oddShape=find(data.Shape<1);
% Vmin=max(data.Volume(oddShape)); %largest volume with impossible shape
% data.case(4).weight=and(data.Shape>2,data.Volume>2*Vmin);
% data.case(4).label='Small and non-elongated volumes filtered out';data.case(4).case='Filtered';
ncase=numel(data.case);
for icase=1:ncase;
    disp(data.case(icase).label);
    data.case(icase).orientation=analyseOrientation(data.XYZ(1,:),...
    data.XYZ(2,:),...
    data.XYZ(3,:),...
    data.case(icase).weight);
end
% 
% 
% % Weighted
% disp(sprintf('Data weighted by volume'));
% data.orientationVolume=analyseOrientation(data.XYZ(1,:),...
%     data.XYZ(2,:),...
%     data.XYZ(3,:),...
%     data.Volume);
% % Weighted by shape
% disp(sprintf('\nData weighted by Shape'));
% data.orientationShape=analyseOrientation(data.XYZ(1,:),...
%     data.XYZ(2,:),...
%     data.XYZ(3,:),...
%     data.Shape);
% % Unweighted
% disp(sprintf('\nUnweighted data'));
% data.orientationNW=analyseOrientation(data.XYZ(1,:),...
%     data.XYZ(2,:),...
%     data.XYZ(3,:),...
%     ones(size(data.Phi)));
% % statsOrientation(data.orientationNW,data.XYZ)

% orientation discreminent plot
for icase=1:ncase;
    OR(1,icase).eigVal=data.case(icase).orientation.eigenValue;
    OR(1,icase).label=sprintf('%s – %s',data.Label,data.case(icase).case);
    switch data.case(icase).case
        case 'Unweighted'
            OR(1,icase).shape='o';OR(1,icase).color=[1,0.647,0]*bright;
        case 'Volume'
            OR(1,icase).shape='s';OR(1,icase).color=[0,0,1]*bright;
        case 'Shape'
            OR(1,icase).shape='d';OR(1,icase).color=[1,0.2,0.6]*bright;
        case 'Filtered'
            OR(1,icase).shape='p';OR(1,icase).color=[0,1,0.3]*bright;
        otherwise
            OR(1,icase).shape='v';OR(1,icase).color=[1,1,1]*bright*0.8;
    end
end  
% plotOrientation(OR,2);

% grid=SchmidtProjection2(data,1,3,30,50);
grid=SchmidtProjectionMT(data,1,3,30,50,OR);
print(1,'-dpdf',sprintf('%s-%s','CQ0705MT',data.Label),'-fillpage');
end

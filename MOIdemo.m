%% Define original object
dl=0.05;lmax=1.5;l=[-lmax:dl:lmax];nl=numel(l);
r=interp1([-1,-0.9,-0.8,-0.7,-0.5,-0.2,0,0.2,0.4,0.6,0.8,1],...
           [0, 0.4, 0.5, 0.2, 0.1, 0.2,0.3,0.3,0.2,0.15,0.1,0],l);
xc=interp1([-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1],...
            [0, 0.2, 0.4, 0.4, 0.2,0,  0,  0,0.2 0.1,0],l);
yc=interp1([-1,-0.8,-0.6,-0.4,-0.2,  0,0.2,0.4,0.6,0.8,1],...
            [0,   0,   0, 0.1, 0.2,0.1,0.3,0.3,0.2,0.1,0.2],l);
[x,y,z]=meshgrid(l);v=x*0;
rho=(x.^2+y.^2).^(1/2);
for il=1:nl;
    rn=r(il);
    % v(:,:,il)=rho(:,:,il)<=r(il);
    v(:,:,il)=((x(:,:,il)-xc(il)).^2+(y(:,:,il)-yc(il)).^2)<=r(il).^2;
end
%%
iv=find(v>0);
pc=[x(iv) y(iv) z(iv)];
% k=convhulln(pc); % find convex hull
k=boundary(pc,0.75); % find convex hull
%
figure(3); clf; 

ht=trisurf(k,pc(:,1),pc(:,2),pc(:,3))
set(ht,'EdgeColor','none','FaceColor','red','FaceAlpha',0.25)
% scatter3(x(iv),y(iv),z(iv));
view(3); axis equal;


% z=(linspace(0,2,100)); x=z*0; y=z*0;
% r=interp1([0,0.5,1,1.5,2]/2,[0,0.3,0.4,0.15,0],z/xmax)
%% Rotate object
yaw=-10; pitch=10;roll=0;vaz=-135;vel=20;
% yaw=30; pitch=50;roll=0;vaz=-25;vel=35;
R=[cosd(pitch),-sind(pitch),0;sind(pitch),cosd(pitch),0;0,0,1]*...
    [cosd(yaw),0,sind(yaw);0,1,0;-sind(yaw),0,cosd(yaw)]*...
    [1,0,0;0,cosd(roll),-sind(roll);0,sind(roll),cosd(roll)];
Z=R*[x(:) y(:) z(:)]';%R*[x(:) y(:) z(:)];

%% Discretize on regular grid
dg=0.05;
xmin=dg*floor(min(Z(1,:))/dg);
xmax=dg*ceil(max(Z(1,:))/dg);
ymin=dg*floor(min(Z(2,:))/dg);
ymax=dg*ceil(max(Z(2,:))/dg);
zmin=dg*floor(min(Z(3,:))/dg);
zmax=dg*ceil(max(Z(3,:))/dg);
[xg,yg,zg]=meshgrid([xmin:dg:xmax],[ymin:dg:ymax],[zmin:dg:zmax]);
vg=griddatan(Z',v(:),[xg(:) yg(:) zg(:)],'nearest');
iv=find(vg>0);
S=[xg(iv) yg(iv) zg(iv)]';

% W=interp3(Z(1,:),Z(2,:),Z(3,:),ones(1,size(Z,2)),xg(:),yg(:),zg(:),'nearest');

%%
nS=size(S,2);
C=mean(S,2);%centroid
P=S-repmat(C,1,nS); %position wrt center;
% R2=sum(diag(P'*P)); %||r_k||^2;
I=sum(diag(P'*P))*eye(3)-P*P'; %Inertia matrix
[V,D]=eig(I);
%scale eigenvector
L=2*(1-diag(D)/sqrt(sum(diag(D).^2)));
VS=V.*repmat(L',3,1);




%%
% pc=[x(iv) y(iv) z(iv)];
% k=convhulln(pc); % find convex hull

K=boundary([P(1,:); P(2,:); P(3,:)]',1); % find convex hull

figure(1); clf; hold on;
% Hs=scatter3(P(1,:),P(2,:),P(3,:));
Ht=trisurf(K,P(1,:),P(2,:),P(3,:))
set(Ht,'EdgeColor','none','FaceColor','red','FaceAlpha',1)
quiver3([0,0,0],[0,0,0],[0,0,0],...^2
    VS(1,:),VS(2,:),VS(3,:),'k','linewidth',2)
quiver3([0,0,0],[0,0,0],[0,0,0],...
    -VS(1,:),-VS(2,:),-VS(3,:),'k','LineWidth',2)


th1=acosd(V(3,1)); la1=atan2d(V(2,1),V(1,1));
na=30;ra=1.2;
lvec=linspace(0,th1,na);
Htheta=plot3(ra*cosd(la1)*sind(lvec),ra*sind(la1)*sind(lvec),ra*cosd(lvec),'LineWidth',4,'color','b');
lvec=linspace(0,la1,na);
Hlambda=plot3(ra*cosd(lvec)*sind(th1),ra*sind(lvec)*sind(th1),lvec*0,'LineWidth',4,'color','b');
HAx=quiver3([0,0],[0,0],[0,0],[0,1.5],[0,0],[1.5,0],'b','LineWidth',1)
Hvert=plot3([0,1,1]*ra*cosd(la1)*sind(th1),[0,1,1]*ra*sind(la1)*sind(th1),[0,0,1]*ra*cosd(th1),'b','LineWidth',1)


view(vaz,vel); axis equal;
set(gca,'visible','off')
Hl=camlight('left');
lighting("gouraud")

dA=30;
% meridians;
A=linspace(0,180,20);
for L=0:dA:(360-dA);
    plot3(cosd(L)*sind(A),sind(L)*sind(A),cosd(A),'k')
end
% parallels;
A=linspace(0,360,40);
for L=[0:dA:(180-dA)];
    plot3(cosd(A)*sind(L),sind(A)*sind(L),cosd(L)*ones(size(A)),'k')
end



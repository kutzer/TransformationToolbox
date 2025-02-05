%% SCRIPT_plotJoint
clear all
close all
clc

%% Create a 3D RRR arm
% Define figure and axes
fig = figure;
axs = axes('Parent',fig,'NextPlot','add','DataAspectRatio',[1 1 1]);
view(axs,3);

% Define scale 
sc = 0.5;

% Define base frame
h_o2a = triad('Parent',axs,'LineWidth',2);

% Joint 1
[h_J1to0,jnt{1}] = plotJoint(h_o2a,'Revolute',[0,0,1],sc);
% Link 1
H_L1toJ1 = Tz(2)*Rx(pi/2);
h_L1toJ1 = triad('Matrix',H_L1toJ1,'Parent',h_J1to0);
lnk1 = plotLink(h_J1to0,H_L1toJ1(1:3,4),(sc/2).*[1.0,0.8]);

% Joint 2
[h_J2toL1,jnt{2}] = plotJoint(h_L1toJ1,'Revolute',[0,0,1],sc);
% Link 2
H_L2toJ2= Tx(3);
h_L2toJ2 = triad('Matrix',H_L2toJ2,'Parent',h_J2toL1);
lnk2 = plotLink(h_J2toL1,H_L2toJ2(1:3,4),(sc/2).*[0.8,0.8]);

% Joint 3
[h_J3toL2,jnt{3}] = plotJoint(h_L2toJ2,'Revolute',[0,0,1],0.5);
% Link 3
H_L3toJ3= Tx(2);
h_L3toJ3 = triad('Matrix',H_L3toJ3,'Parent',h_J3toL2);
lnk3 = plotLink(h_J3toL2,H_L3toJ3(1:3,4),(sc/2).*[0.8,0.0]);

% Show projection of arm
prj_o = plot(h_o2a,nan,nan,'--k','LineWidth',2);
% Highlight end-point of arm
pnt_o = plot(h_o2a,nan,nan,'ok','MarkerFaceColor','k','MarkerSize',8);

% Create update function
jnts = @(q) projectArm(q,jnt,h_L3toJ3,prj_o,pnt_o);

%% Move arm to specific configuration
jnts( deg2rad([20,15,15]) );

%% Clean up plot
hideTriad([h_J1to0,h_L1toJ1,h_J2toL1,h_L2toJ2,h_J3toL2,h_L3toJ3]);
set(axs,'Visible','off','Units','normalized','Position',[0,0,1,1]);
set(fig,'Color',[1,1,1],'Units','Inches','Position',[1,1,9.3,5.4],...
    'PaperSize',[5.4,9.3],'PaperPosition',[0,0,9.3,5.4]);
axis(axs,'tight');

%% Show plane of manipulator
jnts( deg2rad([0,15,15]) );
drawnow
X_lims = [xlim(axs); ylim(axs); zlim(axs)];
X = mean(X_lims,2);
s = norm( diff(X_lims,1,2) )/2;
pln = plotPlane(axs,[0,1,0,0],X,s);
set(pln,'Parent',h_J1to0);

%% Fix axis limits
dX_lims = diff(X_lims,1,2);
dxy_lims = max(dX_lims(1:2));
X_lims(1:2,:) = dxy_lims*repmat([-1,1],2,1);

set(axs,'XLim',X_lims(1,:),'YLim',X_lims(2,:),'ZLim',X_lims(3,:),...
    'XLimMode','manual','YLimMode','manual','ZLimMode','manual');

%% Create videos
vid_3D = VideoWriter('RRR_3D.mp4','MPEG-4');
vid_pZ = VideoWriter('RRR_topdnPln.mp4','MPEG-4');
vid_nY = VideoWriter('RRR_manipPln.mp4','MPEG-4');
open(vid_3D);
open(vid_pZ);
open(vid_nY);

%% Rotate and show views
for phi = linspace(0,360,120)
    % Move arm
    jnts( deg2rad([phi,15,15]) );

    % 3D view
    view(axs,3);
    drawnow
    frm = getframe(fig);
    writeVideo(vid_3D,frm);

    % Negative y-axis view
    H_J1to0 = get(h_J1to0,'Matrix');
    view(axs,-H_J1to0(1:3,2));
    drawnow
    frm = getframe(fig);
    writeVideo(vid_nY,frm);

    % Positive z-axis view
    H_J1to0 = get(h_J1to0,'Matrix');
    view(axs, H_J1to0(1:3,3));
    drawnow
    frm = getframe(fig);
    writeVideo(vid_pZ,frm);
end

close(vid_3D);
close(vid_pZ);
close(vid_nY);

%% Internal functions
% -> Combine joints
function combineJnts(jnt,q)

for i = 1:numel(jnt)
    jnt{i}(q(i));
end

end

% -> Forward kinematics
function H_e2o = fkin(h_L3toJ3)

H_e2o = getAbsoluteTransform(h_L3toJ3);

end

% -> Project arm
function projectArm(q,jnt,h_L3toJ3,prj_o,pnt_o)

combineJnts(jnt,q);

H_e2o = fkin(h_L3toJ3);

X_o(:,1) = zeros(3,1);
X_o(:,2) = [H_e2o(1:2,4); 0];
X_o(:,3) = H_e2o(1:3,4);

set(prj_o,'XData',X_o(1,:),'YData',X_o(2,:),'ZData',X_o(3,:));
set(pnt_o,'XData',H_e2o(1,4),'YData',H_e2o(2,4),'ZData',H_e2o(3,4));

end
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
h_o2a = triad('Parent',axs);

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

% Projections
prj_o = plot(h_o2a,nan,nan,'--k','LineWidth',2);

% Create update function
jnts = @(q) projectArm(q,jnt,h_L3toJ3,prj_o);


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
function projectArm(q,jnt,h_L3toJ3,prj_o)

combineJnts(jnt,q);

H_e2o = fkin(h_L3toJ3);

X_o(:,1) = zeros(3,1);
X_o(:,2) = [H_e2o(1:2,4); 0];
X_o(:,3) = H_e2o(1:3,4);

set(prj_o,'XData',X_o(1,:),'YData',X_o(2,:),'ZData',X_o(3,:));

end
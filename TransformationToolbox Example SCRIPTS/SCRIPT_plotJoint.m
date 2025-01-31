%% SCRIPT_plotJoint
clear all
close all
clc

%% Create a 3D RRR arm
% Define figure and axes
fig = figure;
axs = axes('Parent',fig,'NextPlot','add','DataAspectRatio',[1 1 1]);
view(axs,3);

% Define base frame
h_o2a = triad('Parent',axs);

% Joint 1
[h_J1to0,jnt1] = plotJoint(h_o2a,'Revolute',[0,0,1],0.5);
% Link 1
H_L1toJ1 = Tz(3)*Rx(pi/2);
h_L1toJ1 = triad('Matrix',H_L1toJ1,'Parent',h_J1to0);
lnk1 = plotLink(h_J1to0,H_L1toJ1(1:3,4),[0.5,0.5]);

% Joint 2
[h_J2toL1,jnt2] = plotJoint(h_L1toJ1,'Revolute',[0,0,1],0.5);
% Link 2
H_L2toJ2= Tx(3);
h_L2toJ2 = triad('Matrix',H_L2toJ2,'Parent',h_J2toL1);
lnk2 = plotLink(h_J2toL1,H_L2toJ2(1:3,4),[0.5,0.5]);

% Joint 3
[h_J3toL2,jnt3] = plotJoint(h_L2toJ2,'Revolute',[0,1,0],0.5);
% Link 3
H_L3toJ3= Tx(3);
h_L3toJ3 = triad('Matrix',H_L3toJ3,'Parent',h_J3toL2);
lnk3 = plotLink(h_J3toL2,H_L3toJ3(1:3,4),[0.5,0]);
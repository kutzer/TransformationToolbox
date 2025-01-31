%% SCRIPT_plotJoint
clear all
close all
clc

%% Create a 3D RRR arm
fig = figure;
axs = axes('Parent',fig,'NextPlot','add','DataAspectRatio',[1 1 1]);

h_o2a = triad('Parent',axs);

[h_1to0,jnt1] = plotJoint(h_o2a,'Revolute',[0,0,1],0.5);
h_Lto1 = triad('Matrix',Tz(3),'Parent',h_1to0);
[h_2toL,jnt2] = plotJoint(h_Lto1,'Revolute',[0,1,0],0.5);
h_Lto2 = triad('Matrix',Tx(3),'Parent',h_2toL);
[h_3toL,jnt3] = plotJoint(h_Lto2,'Revolute',[0,1,0],0.5);
h_Lto3 = triad('Matrix',Tx(3),'Parent',h_3toL);
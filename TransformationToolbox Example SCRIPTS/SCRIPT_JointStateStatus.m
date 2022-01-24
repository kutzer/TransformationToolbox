%% SCRIPT_JointStateStatus
clear all
close all
clc

%% Initialize joint configuration
q = 2*pi*rand(6,1) - pi;

%% Create figure for displaying status
fig = figure('Color',[1 1 1]);
axs = axes('Parent',fig,'xlim',[-0.25,0.25],'ylim',[-0.3,0.3],'Visible','off');

msg = sprintf('%9.4f\\\\%9.4f\\\\%9.4f\\\\%9.4f\\\\%9.4f\\\\%9.4f',rad2deg(q));
txt = text(0,0,['$q = \left( \begin{array}{c}',msg,'\end{array} \right)$'],'Interpreter','latex',...
    'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',16);

%% Update text in loop
while true
    % Update q
    q = 2*pi*rand(6,1) - pi;
    
    % Update text in figure
    msg = sprintf('%9.4f\\\\%9.4f\\\\%9.4f\\\\%9.4f\\\\%9.4f\\\\%9.4f',rad2deg(q));
    set(txt,'String',['$q = \left( \begin{array}{c}',msg,'\end{array} \right)$']);
    
    drawnow;
end
    
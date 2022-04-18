%% SCRIPT_approximateFixedKinematics
% Approximate fixed kinematic paramters using collected data
%
%   M. Kutzer, 18Apr2022, USNA

clear all
close all
clc

%% Create data set
L_tru = [151.85; 243.55; 213.20; 131.05; 85.35; 92.10];
q_tru = 2*pi*rand(6,20);

for i = 1:size(q_tru,2)
    H_e2o_tru{i} = fkin(q_tru(:,i),L_tru);
end

%% Define function inputs
%approximateFixedKinematics
fcnH_e2o = @(q,c)fkin(q,c);
H_e2o_ALL = H_e2o_tru;
options.qLims = [-pi,pi].*ones(6,2);
options.q0 = 2*pi*ones(6,1) - pi;
options.cLims = [50,250].*ones(6,2);
options.c0 = 250*rand(6,1);
options.RotationWeight = 20;

%% Approximate fixed kinematics
[c,q,err] = approximateFixedKinematics(fcnH_e2o,H_e2o_ALL,options);
%% SCRIPT_calculateJacobian
% This script tests calculateJacobian.m
%
%   M. Kutzer, 12Apr2022, USNA

clear all
close all
clc

%% Define forward kinematics function
q = sym('q',[6,1]);
c = sym('c',[6,1]);

H_e2o_sym = Rz(q(1))*Tx(c(1))*Rz(q(2))*Ty(c(2))*Tz(q(3))*Rx(c(3))*Ry(q(4))*Tz(c(4))*Rx(q(5))*Ty(c(5))*Rz(q(6))*Tx(c(6));

J = calculateJacobian(q,H_e2o_sym,'Constants',c);
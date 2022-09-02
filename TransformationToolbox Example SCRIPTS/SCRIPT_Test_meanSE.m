%% SCRIPT_Test_meanSE
% Test cases for meanSE.
%
%   M. Kutzer, 02Sep2022, USNA

clear all
close all
clc

%% Define test scale
N    = 10; % Number of values of H
nBad =  5; % Number of bad values of H (if applicable)

%% SE(2) Case 1
% Good values
clc

M = 2;
for i = 1:N
    H{i} = randSE(100,M);
end

muH_c = meanSE(H,'Coupled');
muH_d = meanSE(H,'Decoupled');

%% SE(3) Case 1
% Good values
clc
M = 3;

for i = 1:N
    H{i} = randSE(100,M);
end

muH_c = meanSE(H,'Coupled');
muH_d = meanSE(H,'Decoupled');

%% SE(2) Case 2
% Bad (reflection) values
clc
M = 2;

Reflection = eye(M+1);
Reflection(M,M) = -1;
for i = 1:N
    H{i} = randSE(100,M);
    if i <= nBad
        H{i} = Reflection*H{i};
    end
end

muH_c = meanSE(H,'Coupled');
muH_d = meanSE(H,'Decoupled');

%% SE(3) Case 2
% Bad (reflection) values
clc
M = 3;

Reflection = eye(M+1);
Reflection(M,M) = -1;
for i = 1:N
    H{i} = randSE(100,M);
    if i <= nBad
        H{i} = Reflection*H{i};
    end
end

muH_c = meanSE(H,'Coupled');
muH_d = meanSE(H,'Decoupled');

%% SE(2) Case 3
% Bad (scaled) values
clc
M = 2;

Scale = eye(M+1);
Scale(1:M,1:M) = diag(0.8*ones(1,M));
for i = 1:N
    H{i} = randSE(100,M);
    if i <= nBad
        H{i} = Scale*H{i};
    end
end

muH_c = meanSE(H,'Coupled');
muH_d = meanSE(H,'Decoupled');

%% SE(3) Case 3
% Bad (scaled) values
clc
M = 3;

Scale = eye(M+1);
Scale(1:M,1:M) = diag(0.8*ones(1,M));
for i = 1:N
    H{i} = randSE(100,M);
    if i <= nBad
        H{i} = Scale*H{i};
    end
end

muH_c = meanSE(H,'Coupled');
muH_d = meanSE(H,'Decoupled');

%% Check varargin
muH = meanSE(H,'Coupled',1e-3)
muH = meanSE(H,1e-3,'Coupled')
muH = meanSE(H,false,1e-3,'Coupled')
muH = meanSE(H,1e-3,false,'Coupled')
muH = meanSE(H,1e-3,'Coupled',false)
muH = meanSE(H,'Coupled',true,1e-3)
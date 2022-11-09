function [tfUnique,adj,wAdj] = findUniqueTransforms(H_a2b,dX,dR)
% FINDUNIQUETRANSFORMS finds the set of unique transformations given a
% maximum change in translation and rotation.
%   [tfUnique,wAdj] = findUniqueTransforms(H_a2b,dX,dR)
%
%   Input(s)
%       H_a2b - n-element cell array defining rigid body transformations
%       dX    - [OPTIONAL] positive scalar value defining the minimum
%               tranlation difference to be considered unique. Default
%               value is dX = 5.
%       dR    - [OPTIONAL] positive scalar value defining the minimum 
%               rotation difference to be considered unique. Default value
%               is dR = 0.01.
%
%   Output(s)
%       tfUnique - n-element logical array identifying transformations as
%                  unique (true) or not unique (false)
%       adj      - n x n logical array defining transformations that are
%                  withing dX and dR of one-another
%       wAdj     - n x n x 2 array defining the dX and dR distances between
%                  arrays
%
%   M. Kutzer, 31Oct2022, USNA

%% Check input(s)
narginchk(1,3);

if nargin < 2
    dX = 5;
end

if nargin < 3
    dR = 0.01;
end

% TODO - check inputs

%% Define weighted adjacency
n = numel(H_a2b);

% Define status term
nStatus = 50;
iStatus = unique( round(linspace(1,n,nStatus)) );
msg = 'Calculating weighted adjacency...';
wb = waitbar(0,msg,'Name','findUniqueTransforms.m');
t0_Status = tic;

wAdj = zeros(n,n,2);
for i = 1:n
    H_b2a = invSE(H_a2b{i},true);
    for j = 1:n
        if j < i
            continue
        end
        H_a2a = H_b2a*H_a2b{j};
        
        tErr = norm( H_a2a(1:3,4) );
        rErr = norm( vee( logSO(H_a2a(1:3,1:3),true),true ) );

        wAdj(i,j,1) = tErr;
        wAdj(j,i,1) = tErr;
        wAdj(i,j,2) = rErr;
        wAdj(j,i,2) = rErr;
    end

    if any(i == iStatus)
        prc = i/n;
        ti_Status = toc(t0_Status);
        tf_Status = (1-prc)*ti_Status/prc;
        msg_i = sprintf('~ %d seconds remaining',round(tf_Status));
        waitbar(prc,wb,{msg,msg_i});
        drawnow;
    end
end
delete(wb)

%% Define adjacency
adj = wAdj(:,:,1) < dX & wAdj(:,:,2) < dR;

%% Define tfUnique

% Define status term
msg = 'Calculating unique transforms...';
wb = waitbar(0,msg,'Name','findUniqueTransforms.m');
t0_Status = tic;

tfUnique = false(1,n);
for i = 1:n
    tfChk = adj(i,:);
    tfChk(i) = false; % Remove current node
    if ~any(tfChk & tfUnique)
        tfUnique(i) = true;
    end

    if any(i == iStatus)
        prc = i/n;
        ti_Status = toc(t0_Status);
        tf_Status = (1-prc)*ti_Status/prc;
        msg_i = sprintf('~ %d seconds remaining',round(tf_Status));
        waitbar(prc,wb,{msg,msg_i});
        drawnow;
    end
end
delete(wb);
function SE2 = SE3toSE2(SE3)
% SE3toSE2 converts a 3D homogeneous transformation into a 2D homogeneous
% transformation assuming all rigid motions took place in the xy-plane.
%   SE2 = SE3toSE2(SE3)
%
%   Input(s)
%       SE3 - 4x4 array element of SE(3)
%
%   Output(s)
%       SE2 - 3x3 array element of SE(2)
%
%   NOTES: 
%       (1) This function does not check the validity of the input argument
%           against the properties of SE(3). 
%       (2) This function does not check or confirm that the rotation
%           occurs solely about a z-axis
%
%   M. Kutzer, 26Sep2022, USNA

SE2 = eye(3);
SE2(1:2,1:2) = SE3(1:2,1:2);
SE2(1:2,3)   = SE3(1:2,4);
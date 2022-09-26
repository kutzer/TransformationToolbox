function SE3 = SE2toSE3(SE2)
% SE2toSE3 converts a 2D homogeneous transformation into a 3D homogeneous
% transformation assuming all rigid motions took place in the xy-plane.
%   SE3 = SE2toSE3(SE2)
%
%   Input(s)
%       SE2 - 3x3 array element of SE(2)
%
%   Output(s)
%       SE3 - 4x4 array element of SE(3)
%
%   NOTES: 
%       (1) This function does not check the validity of the input argument
%           against the properties of SE(2). 
%
%   M. Kutzer, 16Feb2016, USNA

% Update(s)
%   26Sep2022 - Updated documentation

SE3 = eye(4);
SE3(1:2,1:2) = SE2(1:2,1:2);
SE3(1:2,4) = SE2(1:2,3);
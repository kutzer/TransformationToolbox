function SE3 = SE2toSE3(SE2)
%SE2toSE3 converts a 2D homogeneous transformation into a 3D homogeneous
%transformation assuming all rigid motions took place in the xy-plane.

SE3(4,4) = 1;
SE3(1:2,1:2) = SE2(1:2,1:2);
SE3(3,3) = 1;
SE3(1:2,4) = SE2(1:2,3);
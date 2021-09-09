function r = logSO(R)
% LOGSO calculates the log of an element of the special orthogonal group.
%
%   M. Kutzer 08Jan2016, USNA

r = logm(R);
if ~isreal(r)
    [Axis,Angle] = SOtoAxisAngle(R);
    v = Axis*Angle;
    r = wedge(v);
end
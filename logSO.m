function r = logSO(R)
% LOGSO calculates the log of an element of the special orthogonal group.

r = logm(R);
if ~isreal(r)
    R3 = R^(1/3);
    r3 = logm(R3);
    r = 3*r3;
end
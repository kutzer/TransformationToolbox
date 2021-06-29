function rbList = rigidBodyListFromIdxList(rbt,idxList)
% RIGIDBODYLISTFROMIDXLIST finds the rigid body list(s) associated with an 
% index list
%   rbList = RIGIDBODYFROMIDXLIST(rbt,IdxList)
%
%   Input(s)
%       rbt     - rigid body tree object
%       idxList - array of indices associated with a unique serial chain of
%                 rigid bodies contained in the rigid body tree. If the
%                 provided index list is incorrect for a given rigid body
%                 tree, the function will throw a warning and return an
%                 empty cell array.
%
%   Output(s)
%       rbList  - cell array containing the serial chain of rigid bodies
%                 associated with the index list
%
%   Examples:
%       rbList = rigidBodyFromIdxList(rbt,[1])
%       rbList{1} = rbt.Base.Children{1};
%
%       rbList = rigidBodyFromIdxList(rb,[2,1])
%       rbList{1} = rbt.Base.Children{2};
%       rbList{2} = rb{1}.Children{1};
%
%       rbList = rigidBodyFromIdxList(rb,[1,3,2])
%       rbList{1} = rbt.Base.Children{1};
%       rbList{2} = rb{1}.Children{3};
%       rbList{3} = rb{2}.Children{2};
%
%   M. Kutzer, 29Jun2021, USNA

%% Check input(s)
narginchk(2,2);
% Check if rbt is a rigidBodyTree object
switch lower(class( rbt ))
    case 'rigidbodytree'
        % Acceptable input
    otherwise
        error('Input must be a valid ribidBodyTree object.');
end
% Check if index list contains only scalars > 0
chk = (idxList == round(idxList)) & (idxList > 0);
if nnz(chk) ~= numel(chk)
    error('Specified index list must contain scalar values greater than 0.');
end

%% Create rigid body tree list
if numel(rbt.Base.Children) >= idxList(1)
    rbList{1} = rbt.Base.Children{ idxList(1) };
else
    warning('RigidBodyTree does not contain an element of idxList(1).');
    rbList = {};
    return
end

for i = 2:numel(idxList)
    if numel(rbList{i-1}.Children) >= idxList(i)
        rbList{i} = rbList{i-1}.Children{ idxList(i) };
    else
        warning('RigidBodyTree does not contain an element of idxList(%d).',i);
        rbList = {};
        return
    end
end
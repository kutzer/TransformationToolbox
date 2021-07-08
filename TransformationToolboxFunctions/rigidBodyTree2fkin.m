function [strs,funcs,links,jlims,jhome] = rigidBodyTree2fkin(rbt)
% RIGIDBODYTREE2FKIN recovers the forward kinematics given a rigid body
% tree object.
%   [str,func,jlims] = RIGIDBODYTREE2FKIN(rbt) recovers forward kinematics
%   given a rigid body tree object. Translations are specified in meters
%   and rotation are specified in radians.
%
%   INPUT(S)
%       rbt - rigid body tree object (see importrobot.m)
%
%   Output(s)
%       strs   - Nx1 cell array, each cell contains a string defining the 
%                forward kinematics using motion primitives of a unique
%                serial chain contained in the rigid body tree object.
%                Motion primitives and functions used:
%               (Rx, Ry, Rz, Tx, Ty, and Tz; see Transformation Toolbox)
%       funcs  - Nx1 cell array, each cell contains a function handle of 
%                forward kinematics given joints defined by an array 
%                q = [q(1); q(2); ...] where q(1) represents the most 
%                proximal joint (closest to the base of the robot). Fixed
%                branches will contain no joint.
%       links  - Nx1 cell array, each cell contains a 1xM cell array with
%                cells that contain string desciptions of the link names
%                for the given serial kinematic chain.
%       jlims  - Nx1 cell array, each cell contains a Kx2 array specifying 
%                the negative and positive joint limits prescribed by 
%                rigid body tree
%       jhome -  Nx1 cell array, each cell contains a Kx1 array specifying 
%                the joint "home" position prescribed by rigid body tree.
%
%   M. Kutzer, 28Jun2021, USNA

debug = false;

%% Check input(s)
narginchk(1,1);
switch lower(class( rbt ))
    case 'rigidbodytree'
        % Acceptable input
    otherwise
        error('Input must be a valid ribidBodyTree object.');
end

%% Expand tree into serial rigid body lists
% Define initial list
kids = rbt.Base.Children;
for i = 1:numel(kids)
    idxLists{i,1} = i;
    rbLists{i,1} = rigidBodyListFromIdxList(rbt,idxLists{i});
    if numel(rbLists{i}{end}.Children) > 0
        fullyExpanded(i,1) = false;
    else
        fullyExpanded(i,1) = true;
    end
end

% Expand all branches
while nnz(fullyExpanded) < numel(fullyExpanded)
    for i = 1:numel(idxLists)
        if ~fullyExpanded(i)
            n = numel(rbLists{i}{end}.Children);
            if n > 0
                idxList = idxLists{i,1};
                idxLists{i,1} = [idxList,1];
                fullyExpanded(i,1) = false;
                rbLists{i} = rigidBodyListFromIdxList(rbt,idxLists{i});
                for j = 2:n
                    idxLists{end+1,1} = [idxList,j];
                    fullyExpanded(end+1,1) = false;
                    rbLists{end+1,1} = rigidBodyListFromIdxList(rbt,idxLists{end});
                end
            else
                fullyExpanded(i,1) = true;
            end
        end
    end
end

% Debug code
if debug
    fprintf('Branch index lists:\n');
    for i = 1:numel(idxLists)
        fprintf('\tBranch %d: [',i);
        n = numel(idxLists{i});
        for j = 1:n
            if j < n
                fprintf('%d, ',idxLists{i}(j));
            else
                fprintf('%d]\n',idxLists{i}(j));
            end
        end
    end
end
            
%% Parse kinematics for each serial rigid body list
% Define zero parameters
rpyZERO = 1e-5;    % define "zero" for angles (radians)
xyzZERO = 1e-5;    % define "zero" for translations (meters)

% Define number of serial rigid body lists
n = numel(rbLists);

% Initialize outputs
strs = cell(n,1);
funcs = cell(n,1);
links = cell(n,1);
jlims = cell(n,1);
jhome = cell(n,1);

% Cycle through each rigid body list
for i = 1:n
    str_i = '';
    links_i = {};
    jlims_i = [];
    jhome_i = [];
    
    rbList = rbLists{i};
    
    % Debug Code
    if debug
        fprintf('Branch %d:\n',i);
    end
    
    % Cycle through each element of the rigid body list
    for j = 1:numel(rbList)
        % Parent transform
        H_j2p = rbList{j}.Joint.JointToParentTransform;
        % Child transform
        H_c2j = rbList{j}.Joint.ChildToJointTransform;
        
        % Recover joint information and limits
        bName = rbList{j}.Name;
        jType = rbList{j}.Joint.Type;
        jAxis = rbList{j}.Joint.JointAxis;
        jLims = rbList{j}.Joint.PositionLimits;
        jHome = rbList{j}.Joint.HomePosition;
        
        % Debug code
        if debug
            fprintf('\tRigid Body %d\n',j);
            fprintf('\t\tName: "%s"\n',bName);
            fprintf('\t\tJoint Type: "%s"\n',jType);
            fprintf('\t\tJoint Axis: [%.5f, %.5f, %.5f]\n',jAxis);
            fprintf('\t\tJoint Limits: [%.5f, %.5f]\n',jLims);
            fprintf('\t\tJoint Home: %.5f\n',jHome);
        end
            
        % Update link name list
        links_i{end+1} = bName;
        
        % Recover "rpy" and "xyz"
        rpy_j2p = rotm2eul(H_j2p(1:3,1:3),'XYZ');
        xyz_j2p = H_j2p(1:3,4).';
        rpy_c2j = rotm2eul(H_c2j(1:3,1:3),'XYZ');
        xyz_c2j = H_c2j(1:3,4).';
        
        % Debug code
        if debug
            fprintf('\t\tJoint relative to Parent:\n')
            fprintf('\t\t\txyz: [%.5f, %.5f, %.5f]\n',xyz_j2p);
            fprintf('\t\t\trpy: [%.5f, %.5f, %.5f]\n',rpy_j2p);
            fprintf('\t\tChild relative to Joint:\n')
            fprintf('\t\t\txyz: [%.5f, %.5f, %.5f]\n',xyz_c2j);
            fprintf('\t\t\trpy: [%.5f, %.5f, %.5f]\n',rpy_c2j);
        end
        
        % Define rigid transformation strings
        % - xyz translation
        str_k{1} = xyz2str(xyz_j2p,xyzZERO);
        str_k{3} = xyz2str(xyz_c2j,xyzZERO);
        % - rpy rotation
        str_k{2} = rpy2str(rpy_j2p,rpyZERO);
        str_k{4} = rpy2str(rpy_c2j,rpyZERO);
        
        % Combine forward kinematics strings for parent
        for k = 1:2
            if ~isempty(str_k{k})
                if isempty(str_i)
                    str_i = sprintf('%s%s',str_i,str_k{k});
                else
                    str_i = sprintf('%s*%s',str_i,str_k{k});
                end
            end
        end
        
        % Add joint transformation
        switch lower( jType )
            case 'fixed'
                % No transformation needs to be added
                str_jnt = '';
            case 'revolute'
                jlims_i(end+1,:) = jLims;
                jhome_i(end+1,:) = jHome;
                % Define revolute joint transformation
                [xyz,sgn] = jointAxis2dir(jAxis,xyzZERO);
                if numel(xyz) == 1
                    str_jnt = sprintf('R%s( %sq(%d) )',xyz,sgn,size(jlims_i,1));
                else
                    str_jnt = sprintf('SOtoSE( Rodrigues(%s,%sq(%d)) )',xyz,sgn,size(jlims_i,1));
                end
            case 'prismatic'
                jlims_i(end+1,:) = jLims;
                jhome_i(end+1,:) = jHome;
                % Define prismatic joint transformation
                [xyz,sgn] = jointAxis2dir(jAxis,xyzZERO);
                if numel(xyz) == 1
                    str_jnt = sprintf('T%s(%sq(%d))',xyz,sgn,size(jlims_i,1));
                else
                    n = significantDigits(jAxis);
                    val_str = ['%.',int2str(n),'f*%s'];
                    evl_str = sprintf('Tx(%s)*Ty(%s)*Tz(%s)',val_str,val_str,val_str);
                    jnt_str = sprintf('q(%d)',size(jlims_i,1));
                    str_jnt = sprintf(evl_str,...
                        jAxis(1),jnt_str,...
                        jAxis(2),jnt_str,...
                        jAxis(3),jnt_str);
                end
            otherwise
                error('Body %d has a joint type "%s" that is not recognized.',i,jType);
        end
        
        % Append joint
        if ~isempty(str_jnt)
            if isempty(str_i)
                str_i = sprintf('%s%s',str_i,str_jnt);
            else
                str_i = sprintf('%s*%s',str_i,str_jnt);
            end
        end
        
        % Combine forward kinematics strings for child
        for k = 3:4
            if ~isempty(str_k{k})
                if isempty(str_i)
                    str_i = sprintf('%s%s',str_i,str_k{k});
                else
                    str_i = sprintf('%s*%s',str,str_k{k});
                end
            end
        end
        
    end % rbList expansion
    
    % Account for empty string
    if isempty(str_i)
        str_i = 'eye(4)';
    end
    
    % Append string, limits, home, and create function
    strs{i,1} = str_i;
    funcs{i,1} = eval( sprintf('@(q)%s',str_i) );
    links{i,1} = links_i;
    jlims{i,1} = jlims_i;
    jhome{i,1} = jhome_i;
end
end

function str = xyz2str(xyz,ZERO)
% XYZ2STR converts an xyz translation into a string of rigid body
% transformation primitives removing transformations associated with ZERO
% translation.
n = significantDigits(ZERO);
str_empt = ['%sT%s(%0.',int2str(n),'f)'];
str_comp = ['%s*T%s(%0.',int2str(n),'f)'];
str = '';
xyz_str = 'xyz';
for i = 1:numel(xyz)
    if abs(xyz(i)) < ZERO
        % Skip transformation
    else
        if isempty(str)
            str = sprintf(str_empt,str,xyz_str(i),xyz(i));
        else
            str = sprintf(str_comp,str,xyz_str(i),xyz(i));
        end
    end
end
end

function str = rpy2str(rpy,ZERO)
% RPY2STR converts a roll/pitch/yaw array into a string of rigid body
% transformation primitives
n = significantDigits(ZERO);
str_empt = ['%sR%s(%0.',int2str(n),'f)'];
str_comp = ['%s*R%s(%0.',int2str(n),'f)'];
str = '';
xyz_str = 'xyz';
for i = 1:numel(rpy)
    c = abs(rpy(i)/pi); % Factor of pi
    s = sign(rpy(i));   % Sign of angle 
    if abs(rpy(i)) < ZERO
        % Skip transformation
    else
        % Check if "c" is a round number or fraction
        fracList = 1:6;
        fracBin = abs( (c.*fracList) - round(c.*fracList) ) < ZERO;
        if nnz(fracBin) > 0
            % Factor of pi
            [num,den] = rat( ...
                round(c.*min(fracList(fracBin))) /...
                min(fracList(fracBin)) );
            if den == 1 && num == 1
                if s > 0
                    str_coef = '';
                else
                    str_coef = '-';
                end
            elseif den == 1
                str_coef = sprintf('%d*',s*c);
            else
                str_coef = sprintf('(%d/%d)*',s*num,den);
            end
            
            if isempty(str)
                str = sprintf('%sR%s( %spi )' ,str,xyz_str(i),str_coef);
            else
                str = sprintf('%s*R%s( %spi )',str,xyz_str(i),str_coef);
            end
        else
            % Decimal value
            if isempty(str)
                str = sprintf(str_empt,str,xyz_str(i),rpy(i));
            else
                str = sprintf(str_comp,str,xyz_str(i),rpy(i));
            end
        end
    end
end
end

function [xyz,sgn] = jointAxis2dir(v,ZERO)
% JOINTAXIS2DIR outputs a string containing x/y/z and sign/scale value
% string associated with the magnitude of rotation, and a unit vector otherwise.
sgn = 1;
nV = norm(v);
if abs( nV - 1 ) > ZERO
    warning('Joint axis is not a unit vector (|v| = %.11f). Appending norm of vector to sign of angle.',nV);
    v = v./nV;
    sgn = sgn * nV;
end

xyz_str = 'xyz';
bin = abs( abs(v) - 1 ) < ZERO;
if nnz(bin) == 1
    xyz = xyz_str(bin);
    sgn = sgn * sign(v(bin));
else
    n = significantDigits(v);
    str_val = ['%0.',int2str(n),'f'];
    nn = numel(v);
    xyz = '[';
    for i = 1:nn
        switch i
            case nn
                str_tmp = sprintf('%s%s]',xyz,str_val);
                xyz = sprintf(str_tmp,v(i));
            otherwise
                str_tmp = sprintf('%s%s, ',xyz,str_val);
                xyz = sprintf(str_tmp,v(i));
        end
    end
end

switch sgn
    case 1
        sgn = '';
    case -1
        sgn = '-';
    otherwise
        n = significantDigits(nV);
        str_val = ['%0.',int2str(n),'f*'];
        sgn = sprintf(str_val,sgn);
end

end

function n = significantDigits(x)
% SIGNIFICANTDIGITS returns the number of significant digits given a value
% or an array of values. For an array of values, the maximum number of
% significan digits is returned.
x = abs(x); %in case of negative numbers
n=0;
for i = 1:numel(x)
    while (floor(x(i)*10^n)~=x(i)*10^n)
        n=n+1;
    end
end
end
%% SCRIPT_DrawCircleWithUR3e
% Draw a circle in the x/z plane using the UR3 manipulator Jacobian
% within a function automatically derived using symbolic variables.
%
%   M. Kutzer, 10Nov2021, USNA

clear all
close all
clc

%% Initialize UR3 simulation (this requires the URSimulationToolbox)
sim = URsim;
sim.Initialize('UR3');

%% Define forward kinematics symbolically
syms q1 q2 q3 q4 q5 q6
q = [q1; q2; q3; q4; q5; q6];
H_e2o_sym = UR_fkin('UR3',q);

%% Derive Jacobian
% funcJ is an anonymous function, you can also just create a function for J
funcJ = calculateJacobian(q,H_e2o_sym);

%% Move arm to random initial configuration
sim.Joints = 2*pi*rand(6,1);

q_now = sim.Joints;                 % Get the current joint configuration
H_e2o_now = UR_fkin('UR3',q_now);   % Get the current end-effector pose

%% Setup plots
plt_Xt(1) = plot(sim.hFrame0,0,0,'m');
plt_Xt(2) = plot(sim.hFrame0,0,0,'dm');

%% Loop through each desired point and move the arm to the desired point
err = 0.01; 
maxStep_xyz = 5;
maxStep_rot = deg2rad(5);
maxStep_q = deg2rad(3);

% Define variable to collect joint position ("inverse kinematic") data
q_all = [];

firstStep = true;
r = 100;
phi_all = linspace(0,2*pi,100); % Discrete values for parameterization
for phi = phi_all
    % Define desired position
    x = r*cos(phi);
    y = 300;
    z = r*sin(phi) + 250;
    % Define desired orientation
    x_hat = [-sin(phi); 0; cos(phi)];
    z_hat = [0; 1; 0];
    y_hat = cross(z_hat,x_hat);
    % Define desired end-effector pose
    H_e2o_des = [x_hat, y_hat, z_hat, [x;y;z]; 0,0,0,1];
    
    h = triad('Parent',sim.hFrame0,'Matrix',H_e2o_des,'Scale',30);
    
    % Move the arm to the desired configuration
    while true
        delta_H_e2o = invSE(H_e2o_now) * H_e2o_des;

        % Calculate delta_X
        delta_X(1:3,1) = delta_H_e2o(1:3,4);
        delta_X(4:6,1) = vee( logm(delta_H_e2o(1:3,1:3)),'fast');
        
        % Reference delta_X to the base frame
        R_e2o_now = H_e2o_now(1:3,1:3);
        delta_X(1:3,1) = R_e2o_now*delta_X(1:3,1);
        %delta_X(4:6,1) = R_e2o_now*delta_X(4:6,1);
        
        % Plot delta_X
        Xo = H_e2o_now(1:3,4);
        Xf = Xo + delta_X(1:3);
        set(plt_Xt(1),'XData',[Xo(1),Xf(1)],'YData',[Xo(2),Xf(2)],'ZData',[Xo(3),Xf(3)]);
        set(plt_Xt(2),'XData',Xf(1),'YData',Xf(2),'ZData',Xf(3));
        
        % Bound the xyz-translational velocity
        if norm(delta_X(1:3)) > maxStep_xyz
            delta_X(1:3) = maxStep_xyz*( delta_X(1:3)./norm(delta_X(1:3)) );
        end
        % Bound the rotational velocity
        if norm(delta_X(4:6)) > maxStep_rot
            delta_X(4:6) = maxStep_rot*( delta_X(4:6)./norm(delta_X(4:6)) );
        end
        
        % Calculate delta_q
        J = funcJ(q_now);
        %delta_q = ( J^(-1) )*delta_X;
        delta_q = pinv( J )*delta_X;
        
        % Limit change in joint configuration
        if max( abs(delta_q) ) > maxStep_q
            delta_q = maxStep_q*( delta_q ./ max( abs(delta_q) ) );
        end
        fprintf('[')
        fprintf('%.4f ',delta_q);
        fprintf(']\n');
        
        % Calculate updated q
        q_now = q_now + delta_q;
        % Calculate updated H_e2o
        H_e2o_now = UR_fkin('UR3',q_now);
        
        % Update visualization
        sim.Joints = q_now;
        drawnow;
        
        % Check if we are sufficiently close to the desired waypoint
        X_e2o_now = H_e2o_now(1:3,4); 
        X_e2o_des = H_e2o_des(1:3,4);
        if norm( X_e2o_now - X_e2o_des ) < err
            fprintf('Waypoint phi = %0.4f achieved\n', phi)
            break
        end
    end
    
    % Collect joint configuration data
    q_all(:,end+1) = q_now;
end 

%% Plot data
figure;
axes;
hold on
xlabel('Iteration');
ylabel('Joint Angle (deg)');
for i = 1:size(q_all,1)
    plot( rad2deg(q_all(i,:)) );
end
legend('\theta_1','\theta_2','\theta_3','\theta_4','\theta_5','\theta_6');
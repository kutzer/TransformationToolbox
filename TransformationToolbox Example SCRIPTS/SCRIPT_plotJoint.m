%% SCRIPT_plotJoint
clear all
close all
clc

%% Create a 3D RRR arm
% Define figure and axes
fig = figure;
axs = axes('Parent',fig,'NextPlot','add','DataAspectRatio',[1 1 1]);
view(axs,3);

% Define scale
sc = 0.5;

% Define base frame
h_o2a = triad('Parent',axs,'LineWidth',2);

% Joint 1
[h_J1to0,jnt{1}] = plotJoint(h_o2a,'Revolute',[0,0,1],sc);
% Link 1
H_L1toJ1 = Tz(2)*Rx(pi/2);
h_L1toJ1 = triad('Matrix',H_L1toJ1,'Parent',h_J1to0);
lnk1 = plotLink(h_J1to0,H_L1toJ1(1:3,4),(sc/2).*[1.0,0.8]);

% Joint 2
[h_J2toL1,jnt{2}] = plotJoint(h_L1toJ1,'Revolute',[0,0,1],sc);
% Link 2
H_L2toJ2= Tx(3);
h_L2toJ2 = triad('Matrix',H_L2toJ2,'Parent',h_J2toL1);
lnk2 = plotLink(h_J2toL1,H_L2toJ2(1:3,4),(sc/2).*[0.8,0.8]);

% Joint 3
[h_J3toL2,jnt{3}] = plotJoint(h_L2toJ2,'Revolute',[0,0,1],0.5);
% Link 3
H_L3toJ3= Tx(2);
h_L3toJ3 = triad('Matrix',H_L3toJ3,'Parent',h_J3toL2);
lnk3 = plotLink(h_J3toL2,H_L3toJ3(1:3,4),(sc/2).*[0.8,0.0]);

% Show projection of arm
prj_o = plot(h_o2a,nan,nan,'--k','LineWidth',2);
% Highlight end-point of arm
pnt_o = plot(h_o2a,nan,nan,'ok','MarkerFaceColor','k','MarkerSize',8);

% Create update function
jnts = @(q) projectArm(q,jnt,h_L3toJ3,prj_o,pnt_o);

%% Move arm to specific configuration
jnts( deg2rad([20,15,15]) );

%% Clean up plot
hideTriad([h_J1to0,h_L1toJ1,h_J2toL1,h_L2toJ2,h_J3toL2,h_L3toJ3]);
set(axs,'Visible','off','Units','normalized','Position',[0,0,1,1]);
set(fig,'Color',[1,1,1],'Units','Inches','Position',[1,1,9.3,5.4],...
    'PaperSize',[5.4,9.3],'PaperPosition',[0,0,9.3,5.4]);
axis(axs,'tight');

%% Show plane of manipulator
jnts( deg2rad([0,15,15]) );
drawnow
X_lims = [xlim(axs); ylim(axs); zlim(axs)];
X = mean(X_lims,2);
s = norm( diff(X_lims,1,2) )/2;
pln = plotPlane(axs,[0,1,0,0],X,s);
set(pln,'Parent',h_J1to0);

%% Fix axis limits
dX_lims = diff(X_lims,1,2);
dxy_lims = max(dX_lims(1:2));
X_lims(1:2,:) = dxy_lims*repmat([-1,1],2,1);

set(axs,'XLim',X_lims(1,:),'YLim',X_lims(2,:),'ZLim',X_lims(3,:),...
    'XLimMode','manual','YLimMode','manual','ZLimMode','manual');

%% Create videos
fname_3D = 'RRR_3D';
fname_pZ = 'RRR_topdnPln';
fname_nY = 'RRR_manipPln';
vid_3D = VideoWriter([fname_3D,'.mp4'],'MPEG-4');
vid_pZ = VideoWriter([fname_pZ,'.mp4'],'MPEG-4');
vid_nY = VideoWriter([fname_nY,'.mp4'],'MPEG-4');
open(vid_3D);
open(vid_pZ);
open(vid_nY);

%% Rotate and show views
xx = {};
x_m_star = [];
for phi = linspace(0,360,120)
    % Reset axes position
    set(axs,'Position',[0,0,1,1]);

    % Move arm
    jnts( deg2rad([phi,15,15]) );

    % 3D view
    view(axs,3);
    if numel(xx) < 1
        xx{1}(1,:) = xlim(axs);
        xx{1}(2,:) = ylim(axs);
        xx{1}(3,:) = zlim(axs);
    else
        xlim(axs,xx{1}(1,:));
        ylim(axs,xx{1}(2,:));
        zlim(axs,xx{1}(3,:));
    end
    drawnow
    frm = getframe(fig);
    writeVideo(vid_3D,frm);

    % Positive z-axis view (Top-Down)
    H_J1to0 = get(h_J1to0,'Matrix');
    view(axs, H_J1to0(1:3,3));
    if numel(xx) < 2
        xx{2}(1,:) = xlim(axs);
        xx{2}(2,:) = ylim(axs);
        xx{2}(3,:) = zlim(axs);
    else
        xlim(axs,xx{2}(1,:));
        ylim(axs,xx{2}(2,:));
        zlim(axs,xx{2}(3,:));
    end
    drawnow
    frm = getframe(fig);
    writeVideo(vid_pZ,frm);

    % Negative y-axis view (Plane of the Manipulator
    mag = 1;
    while true
        % Update axes position
        set(axs,'Position',[((1-mag)/2)*ones(1,2), mag*ones(1,2)]);

        H_J1to0 = get(h_J1to0,'Matrix');
        view(axs,-H_J1to0(1:3,2));
        if numel(xx) < 3
            xx{3}(1,:) = xlim(axs);
            xx{3}(2,:) = ylim(axs);
            xx{3}(3,:) = zlim(axs);
        else
            xlim(axs,xx{3}(1,:));
            ylim(axs,xx{3}(2,:));
            zlim(axs,xx{3}(3,:));
        end
        %axis(axs,'tight');
        drawnow
        frm = getframe(fig);

        % Process image
        im = frm.cdata;
        im = rgb2gray(im);

        [~,x_m] = find(im==0);
        x_m = max(x_m);

        if isempty(x_m_star)
            x_m_star = x_m;
            break
        end

        if x_m == x_m_star
            fprintf('%.3f\n',mag);
            break
        else
            if mag == 1
                fprintf('%03d, %03d: ',x_m_star,x_m);
                mag = x_m_star/x_m;
            end
            %mag = x_m/x_m_star;
            if x_m > x_m_star
                mag = mag - 0.001;
            else
                mag = mag + 0.001;
            end
        end
    end
    writeVideo(vid_nY,frm);

end

close(vid_3D);
close(vid_pZ);
close(vid_nY);

%% Take snapshots
phi = 30;

% Reset axes position
set(axs,'Position',[0,0,1,1]);

% Move arm
jnts( deg2rad([phi,15,15]) );

% 3D view
view(axs,3);
    xlim(axs,xx{1}(1,:));
    ylim(axs,xx{1}(2,:));
    zlim(axs,xx{1}(3,:));
drawnow
frm = getframe(fig);
imwrite(frm.cdata,[fname_3D,'.png'],'png');

% Positive z-axis view (Top-Down)
H_J1to0 = get(h_J1to0,'Matrix');
view(axs, H_J1to0(1:3,3));
if numel(xx) < 2
    xx{2}(1,:) = xlim(axs);
    xx{2}(2,:) = ylim(axs);
    xx{2}(3,:) = zlim(axs);
else
    xlim(axs,xx{2}(1,:));
    ylim(axs,xx{2}(2,:));
    zlim(axs,xx{2}(3,:));
end
drawnow
frm = getframe(fig);
imwrite(frm.cdata,[fname_pZ,'.png'],'png');

% Negative y-axis view (Plane of the Manipulator
mag = 1;
while true
    % Update axes position
    set(axs,'Position',[((1-mag)/2)*ones(1,2), mag*ones(1,2)]);

    H_J1to0 = get(h_J1to0,'Matrix');
    view(axs,-H_J1to0(1:3,2));
    if numel(xx) < 3
        xx{3}(1,:) = xlim(axs);
        xx{3}(2,:) = ylim(axs);
        xx{3}(3,:) = zlim(axs);
    else
        xlim(axs,xx{3}(1,:));
        ylim(axs,xx{3}(2,:));
        zlim(axs,xx{3}(3,:));
    end
    %axis(axs,'tight');
    drawnow
    frm = getframe(fig);

    % Process image
    im = frm.cdata;
    im = rgb2gray(im);

    [~,x_m] = find(im==0);
    x_m = max(x_m);

    if isempty(x_m_star)
        x_m_star = x_m;
        break
    end

    if x_m == x_m_star
        fprintf('%.3f\n',mag);
        break
    else
        if mag == 1
            fprintf('%03d, %03d: ',x_m_star,x_m);
            mag = x_m_star/x_m;
        end
        %mag = x_m/x_m_star;
        if x_m > x_m_star
            mag = mag - 0.001;
        else
            mag = mag + 0.001;
        end
    end
end
imwrite(frm.cdata,[fname_nY,'.png'],'png');

%% Internal functions
% -> Combine joints
function combineJnts(jnt,q)

for i = 1:numel(jnt)
    jnt{i}(q(i));
end

end

% -> Forward kinematics
function H_e2o = fkin(h_L3toJ3)

H_e2o = getAbsoluteTransform(h_L3toJ3);

end

% -> Project arm
function projectArm(q,jnt,h_L3toJ3,prj_o,pnt_o)

combineJnts(jnt,q);

H_e2o = fkin(h_L3toJ3);

X_o(:,1) = zeros(3,1);
X_o(:,2) = [H_e2o(1:2,4); 0];
X_o(:,3) = H_e2o(1:3,4);

set(prj_o,'XData',X_o(1,:),'YData',X_o(2,:),'ZData',X_o(3,:));
set(pnt_o,'XData',H_e2o(1,4),'YData',H_e2o(2,4),'ZData',H_e2o(3,4));

end
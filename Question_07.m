%% QUESTION 06
clear all

%% CONSTANTS
global g;       % source strengths indexed by plate clockwise
global pt;      % list of (x,y) cooreds for points that make up the body
global Uinf;    % horizontal component of ambient flow velocity
global Vinf;    % verticle component of ambient flow velocity
global Xmj;     % X coords for mid point of planels
global Ymj;     % Y coords for mid point of planels
global S;       % length of panels indexed by plate clockwise
global a;       % in {local} start of panel
global b;       % in {local} end of panel


%% BUILD POINTS FOR PANELS
R = 1;            % radius
N = 64;           % number of point
c = 0.95;         % jowkow circle raius
xs = -0.04875;    % x offset
ys = 0.05;        % y offset
aoa = deg2rad(25); % angle of attack
offset = mod(atan2(-ys,(abs(xs)+c)),2*pi);
angles = flip(linspace(offset,2*pi+offset,N+1));
pt_og =R*[cos(angles);sin(angles)];

% zowkowsky time, aerfoils in disguise as circles
pt_z = complex(pt_og(1,:),pt_og(2,:));
xs_z = xs + 0*1i;
ys_z = 0 + ys*1i;

% trasformation 01
pt_zcs = pt_z + xs_z + ys_z; % offset for camber
% trasformation 02
pt_zj = pt_zcs + c^2./pt_zcs; % airfoil

pt(1,:) = real(pt_zj);
pt(2,:) = imag(pt_zj);

%% ESTABLISH AXIS HANDLES
figure('Name','Streamlines','NumberTitle','off');
stream_lines = axes;
hold(stream_lines,'on');
axis(stream_lines,'equal');
axis(stream_lines,[-5 5 -3 3]);
caxis(stream_lines,[0.5 1.5]);
colormap(stream_lines,'jet');
xlabel(stream_lines,'x');
ylabel(stream_lines,'y');
shading(stream_lines,'flat');
title1 = sprintf('Streamlines and Velocity gradient for a cambered airfoil in a cross flow \napproximated with %i panels',N);
title(stream_lines,title1)



%% SET UP CONSTANTS
Uinf_actual = 1;   % background flow velocity
Vinf_actual = 0;
% rotate flow instead of airfoil
Uinf = Uinf_actual*cos(aoa);
Vinf = Uinf_actual*sin(aoa);

for plate = 1 : N
    Xmj(plate) = (pt(1,plate+1) - pt(1,plate))/2 + pt(1,plate); % centre points
    Ymj(plate) = (pt(2,plate+1) - pt(2,plate))/2 + pt(2,plate); % centre points
    S(plate) = sqrt( (pt(1,plate+1) - pt(1,plate))^2 + (pt(2,plate+1) - pt(2,plate))^2 ); % length of each panel
end

a = -S./2;
b = S./2;

%% GET THE COEFICIENTS FOR NORMAL VELOCITY INDUCED BY PLATE
psi = zeros(N+1,N+1);
psi_inf_local = zeros(1,N+1);

for ii = 1 : N
    for jj = 1 : N
        Xi = [pt(1,ii) pt(1,ii+1)];
        Yi = [pt(2,ii) pt(2,ii+1)];
        Xj = [pt(1,jj) pt(1,jj+1)];
        Yj = [pt(2,jj) pt(2,jj+1)];
        [psi(ii,jj), psi_inf_local(ii)] = GET_vortex_induced_norm_coeff(Xi,Yi,Xj,Yj);
        % psi the totale sum of flow induced at the centre of each plate
        % psi_inf_local the normal component of ambient flow at the centre
        % of the plate
    end
end

% satisfy the cutter condition
psi(:,N+1)=1;
psi(N+1,1)=1;
psi(N+1,N)=1;
psi(N+1,2:(N-1))=0;
psi(N+1,N+1)=0;
psi_inf_local(1,N+1)=0;
g = (psi\(psi_inf_local'))'; % over constrained circulation
g(size(g,2)) = []; % delete the last element C

%% BUILD ST-REAMLINES
dt = .05;
t0 = 0;
tf = 25;
tick = 1;

end_y = 1;
start_y = -8;
step_y = 0.2;
xp = -6;

for yp = start_y : step_y : end_y
    XY = RK_4_2vars_vortex(dt,t0,tf,xp,yp);
    x_pos(tick,:) = XY(1,:); % x_positions along stream lines
    y_pos(tick,:) = XY(2,:); % y positions along stream lines
    dx(tick,:) = XY(3,:);    % corresponding velocity in the dx/dt = u direction
    dy(tick,:) = XY(4,:);    % and dy/dt = v direction
    tick = tick + 1;
end


% add streamline from trailing edge
XY_trail = RK_4_2vars_vortex(dt,t0,tf,pt(1,1)+.001,pt(2,1));


% ROTATE BACK TO HORIZONTAL
pt_z_actual = pt_zj.*exp(1)^(-1i*aoa); % rotation
z_coords = x_pos + y_pos*1i;
z_dxdy = dx + dy*1i;
z_trail = XY_trail(1,:) + XY_trail(2,:)*1i;

z_dxdy = z_dxdy.*exp(1)^(-1i*aoa); % rotation
z_coords = z_coords.*exp(1)^(-1i*aoa); % rotation
z_trail = z_trail.*exp(1)^(-1i*aoa); % rotation

pt_actual(1,:) = real(pt_z_actual); % rotation
pt_actual(2,:) = imag(pt_z_actual); % rotation
x_pos = real(z_coords);
y_pos = imag(z_coords);
dx = real(z_dxdy);
dy = imag(z_dxdy);
XY_trail(1,:) = real(z_trail);
XY_trail(2,:) = imag(z_trail);



%% PLOTS
[hC hC] = contourf(stream_lines,x_pos,y_pos,sqrt(dx.^2+dy.^2),100);
%** note velocity info is taken directly from streamlines **%
set(hC,'LineStyle','none'); % turn on contour boundry lines

for line = 1 : size(x_pos,1)
    plot(stream_lines,x_pos(line,:),y_pos(line,:),'k')
end
fill(stream_lines, pt_actual(1,:), pt_actual(2,:),'k'); % plot the airfoil
plot(stream_lines,XY_trail(1,:),XY_trail(2,:),'m') % plot the streamline @ trailing edge


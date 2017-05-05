%% QUESTION 01
clear all

%% CONSTANTS
global q;       % source strengths indexed by plate clockwise
global pt;      % list of (x,y) cooreds for points that make up the body
global Uinf;    % horizontal component of ambient flow velocity
global Vinf;    % verticle component of ambient flow velocity
global Xmj;     % X coords for mid point of planels
global Ymj;     % Y coords for mid point of planels
global S;       % length of panels indexed by plate clockwise
global a;       % in {local} start of panel
global b;       % in {local} end of panel


%% BUILD POINTS FOR PANELS
R = 1;  %radius
N = 64;  % number of point
angles = flip(linspace(pi/64,2*pi+pi/64,N+1));
pt =R*[cos(angles);sin(angles)];
%pt = flip(pt);


%% ESTABLISH AXIS HANDLES
figure('Name','Streamlines','NumberTitle','off');
stream_lines = axes;
hold(stream_lines,'on');
axis(stream_lines,'equal')
xlabel(stream_lines,'x')
ylabel(stream_lines,'y')
title1 = sprintf('Streamlines and Velocity gradient for an cylinder in a cross flow \napproximated with %i panels',N); 
title(stream_lines,title1)  
daspect(stream_lines,[1 1 1])
colorbar(stream_lines);

figure('Name','Pressure Distribution','NumberTitle','off');
pressure = axes;
hold(pressure,'on');
axis(pressure,'equal')
xlabel(pressure,'theta')
ylabel(pressure,'C_p')
title1 = sprintf('Pressure Coeficient distribution'); 
title(pressure,title1)
axis(pressure,[pi/64 2*pi -3.5 1])

%% SET UP CONSTANTS
Uinf = 1;   % background flow velocity 
Vinf = 0;

for z = 1 : N
    Xmj(z) = (pt(1,z+1) - pt(1,z))/2 + pt(1,z); % centre points
    Ymj(z) = (pt(2,z+1) - pt(2,z))/2 + pt(2,z); % centre points
    S(z) = sqrt( (pt(1,z+1) - pt(1,z))^2 + (pt(2,z+1) - pt(2,z))^2 ); % length of each panel
end

a = -S./2;
b = S./2;

%% GET THE COEFICIENTS FOR NORMAL VELOCITY INDUCED BY PLATE
v = zeros(N,N);
u = zeros(N,N);
V = zeros(1,N);
U = zeros(1,N);

for i = 1 : N
    for j = 1 : N        
        % hard code trigger to catch plates self inducing velocity
        if i == j
            v(i,j) = 0.5;
        else
            Xi = [pt(1,i) pt(1,i+1)];
            Yi = [pt(2,i) pt(2,i+1)];
            Xj = [pt(1,j) pt(1,j+1)];
            Yj = [pt(2,j) pt(2,j+1)];
            [v(i,j),u(i,j), V(i), U(i)] = GET_induced_norm_coeff(Xi,Yi,Xj,Yj,Uinf);            
            % v(i,j) the normal velocity induced by plate j at plate i
            % u(i,j) the tangential velocity induced by plate j at plate i
            % V the normal velcoity from the background flow at plate i
            % U the tangential component of the background flow at plat i
        end       
    end
end

q = (v\(V)')'; % 1 by 8 vector of plate sources strengths to create solid boundary

%% BUILD STREAMLINES 
dt = 0.1;
t0 = 0;
tf = 7;
start_x = -3;

% get the seperatrix seperatly
sep_rhs = RK_4_2vars(dt,t0,tf,0.93,0);
sep_lhs = RK_4_2vars(dt,t0,tf,-3,0);
sep_rhs_short = sep_rhs(:,1:2:69); % take every second value
sep_lhs_short = sep_lhs(:,1:2:end);

% get the other streamlines
tick = 1;
for start_y = -2.25 : 0.25 : 2.25
    
    if start_y == 0
    % messy but this concatinates the seperatrix on either side of the geometry    
    x_pos(tick,:) = [sep_lhs_short(1,:),sep_rhs_short(1,:)];
    y_pos(tick,:) = [sep_lhs_short(2,:),sep_rhs_short(2,:)];
    dx(tick,:) =  [sep_lhs_short(3,:),sep_rhs_short(3,:)];
    dy(tick,:) =  [sep_lhs_short(4,:),sep_rhs_short(4,:)];
            
    tick = tick + 1;
    else
    
    XY = RK_4_2vars(dt,t0,tf,start_x,start_y);
    %plot(stream_lines,XY(1,:),XY(2,:),'g')
    
    x_pos(tick,:) = XY(1,:); % x_positions along stream lines
    y_pos(tick,:) = XY(2,:); % y positions along stream lines
    dx(tick,:) = XY(3,:);    % corresponding velocity in the dx/dt = u direction
    dy(tick,:) = XY(4,:);    % and dy/dt = v direction
    tick = tick + 1;
    end
end



%% PLOTS
[hC hC] = contourf(stream_lines,x_pos,y_pos,sqrt(dx.^2+dy.^2),100); 
%** note velocity info is taken directly from streamlines **%
set(hC,'LineStyle','none');
shading(stream_lines,'flat')
colormap(stream_lines,'spring'); % make it different to the fig in the notes...

for line = 1 : size(x_pos,1)
    plot(stream_lines,x_pos(line,:),y_pos(line,:),'k')
end
fill(stream_lines, pt(1,:),pt(2,:),'c'); % really different...
axis(stream_lines,[-3 3 -2.2 2.2])


%% END QUESTION 01

%% QUESTION 02
% we have u_i the coefficent in the normal direction, q the source strength 
% and Uinf the normal compoent due to the ambient flow
Uc_i = q*u + U;
Cp_i = 1 - (Uc_i/Uinf).^2;
t2 = [pi/64:0.001:2*pi+pi/64];
Cp_act = 1-(4*(sin(t2)).^2);
plot(pressure,t2,Cp_act,angles(1:end-1)-pi/64,Cp_i, 'r:o') 
% note pi/8 offset is required from where panel i=1 is located

%% END QUESTION 02


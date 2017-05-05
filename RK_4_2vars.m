function [ XYT ] = RK_4_2vars(dt,t0,tf,xi,yi)
% 4TH DEGREE RUNGE CUTTER 
% 2 VARIABLES %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUT
%   - dt, time step
%   - t0, start
%   - tf, end
%   - xi, initial x
%   - yi, initial y

% OUTPUT
%   - XY = [x coords] 
%          [y coords]
%          [u at correspoinding point (x,y)]
%          [v at correspoinding point (x,y)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIAL CONDITIONS
t = t0:dt:tf;
x = zeros(size(t));
y = zeros(size(t));
x(1) = xi;
y(1) = yi;
dx = zeros(size(t));
dy = zeros(size(t));

% EVELUAITON
for i = 1:1:length(t)-1   
          
    [k11, k12] = dxdt_dydt( x(i),            y(i)            );
    [k21, k22] = dxdt_dydt( x(i) + k11*dt/2, y(i) + k12*dt/2 );
    [k31, k32] = dxdt_dydt( x(i) + k21*dt/2, y(i) + k22*dt/2 );
    [k41, k42] = dxdt_dydt( x(i) + k31*dt,   y(i) + k32*dt   );
                
    % x,y step     
    dx(i) = (1/6*(k11+k41) + 1/3*(k21+k31));
    dy(i) = (1/6*(k12+k42) + 1/3*(k22+k32));
    
    % next point
    x(i+1) = x(i) + dx(i)*dt;
    y(i+1) = y(i) + dy(i)*dt;
end

XYT = [x;y;dx;dy]; %[x;y,;t];

%% end RK_4
end

%% DEFINE YOUR EQUATIONS HERE:
function [dxdt, dydt] = dxdt_dydt(xp,yp)

    % GLOBAL VARIABLES
    global q;
    global pt;
    global Uinf;
    global Vinf;
    global Xmj;
    global Ymj;
    global a;
    global b;
    
    % look at evert panel
    for z = 1 : size(pt,2)-1   
        % angles local to this point (xp. yp) and panel # z
        beta(z) = atan2( (yp - Ymj(z)) , (xp - Xmj(z)) );
        phi(z) = atan2( (pt(2,z+1) - pt(2,z)) , (pt(1,z+1) - pt(1,z)) );
        omega(z) = beta(z) - phi(z);        
    end
   
    
    r = sqrt( (xp - Xmj).^2 + (yp - Ymj).^2 );
 
    x_o_prime = r.*cos(omega);
    y_o_prime = r.*sin(omega);
    
    %% need to sum over every plate... 
    u_prime = q./(2*pi).*(...
        (0.5.*log(y_o_prime.^2 + a.^2 -2.*x_o_prime.*a + x_o_prime.^2)) -...
        (0.5.*log(y_o_prime.^2 + b.^2 -2.*x_o_prime.*b + x_o_prime.^2)));   
    
    v_prime = q./(2.*pi).*( atan((b-x_o_prime)./y_o_prime) -...
                            atan((a-x_o_prime)./y_o_prime));                                        
                           
    dxdt = u_prime.*cos(phi) - v_prime.*sin(phi);
    dydt = u_prime.*sin(phi) + v_prime.*cos(phi);    
    dxdt = sum(dxdt) + Uinf;
    dydt = sum(dydt) + Vinf;
    %% END dxdt
end

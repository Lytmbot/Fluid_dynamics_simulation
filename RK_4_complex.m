function [ z ] = RK_4_complex(dt,t0,tf,xi,yi)
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
%   - z = x(t) + iy(t) values in 1D matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIAL CONDITIONS
t = t0:dt:tf;
z = zeros(size(t));
z(1) = xi + 1i*yi;

% EVELUAITON
for j = 1:1:length(t)-1   

    k1 = dwdz(z(j));        
    k2 = dwdz(z(j) + k1*dt/2);    
    k3 = dwdz(z(j) + k2*dt/2);    
    k4 = dwdz(z(j) + k3*dt);      
    k1 = conj(k1);
    k2 = conj(k2);
    k3 = conj(k3);
    k4 = conj(k4);

    z(j+1) = z(j) + ( (1/6)*(k1+k4) + (1/3)*(k2+k3) )*dt;        

end

%% end RK_4
end

%% DEFINE YOUR EQUATIONS HERE:
function [dwdz] = dwdz(z)
    global a;
    global Uinf;
    global TAO;

    dwdz = Uinf*(1-(a^2/z^2)) + 1i*TAO/(2*pi*z);
    %% END dwdz
end

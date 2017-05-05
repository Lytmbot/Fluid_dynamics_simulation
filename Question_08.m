%% BUILD PLOTS
% clear all
% close all

%% INITIAL CONDITIONS
global a;
global Uinf;
global aoa;
global yo;
global xo;
global c;
global TAO;



aoa = -deg2rad(7);
xo = -0.04875;
yo = 0.05*1i;
c  = 0.95;
a  = 1;
Uinf = 1;
TAO = 4*pi*Uinf*a*sin(asin(abs(yo)/a) - aoa);



%% MAKE A CYLINDER TO TRANSFORM
theta = 0:0.1:2*pi+.1;
radius = 1;
z_circle = a.*cos(theta) + 1i.*a.*sin(theta);


%% STREAMLINES COMPLETE
dt = .05;
t0 = 0;
tf = 15;

end_y = 3;
start_y = -4;
step_y = 0.2;
xp = -5.5;

z  = ones(size(xs,2),tf/dt+1);
z1 = ones(size(xs,2),tf/dt+1);
z2 = ones(size(xs,2),tf/dt+1);
z3 = ones(size(xs,2),tf/dt+1);
z4 = ones(size(xs,2),tf/dt+1);
size(z4)

%% EVELUATION
tick = 1;
for yp = start_y : step_y : end_y
    z(tick,:) = RK_4_complex(dt,t0,tf,xp,yp);
    z1(tick,:) = exp(-1i*aoa).*z(tick,:);
    z2(tick,:) = z1(tick,:) + xo + yo;
    z3(tick,:) = z2(tick,:) + c^2./z2(tick,:);
    z4(tick,:) = exp(1i*aoa).*z3(tick,:);
    plot(stream_lines, real(z4(tick,:)),imag(z4(tick,:)),'b--')
    tick = tick + 1;
end



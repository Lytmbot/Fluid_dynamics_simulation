function [psi, psi_inf_local] = GET_vortex_induced_norm_coeff( Xi, Yi, Xj, Yj )
%% DESRIPTION
%  takes the location of two plates i and j and returns the coffecicent
%  for the induced flow normal to the centre of plate i from j

%% IN
%       start and end ploints for two plates
%       Uing ambient flow velocity

%% OUT
%       psi
%
%       psi_inf_local

global Uinf;
global Vinf;

%set source strength to unity to isolate coefficicent
g_j = 1;

% panel i
Xmi=0.5*(Xi(2)+Xi(1)); % midpoint
Ymi=0.5*(Yi(2)+Yi(1));

% panel j
Xmj=0.5*(Xj(2)+Xj(1)); % midpoint
Ymj=0.5*(Yj(2)+Yj(1));
Phi_j=atan2((Yj(2)-Yj(1)),(Xj(2)-Xj(1))); %plate j angle (eqn 23)

rij = sqrt((Xmj - Xmi).^2 + (Ymj - Ymi).^2); % (eqn 22)

beta = atan2((Ymi - Ymj),(Xmi - Xmj)); % (eqn 25)
omega = beta - Phi_j; % (eqn 26)

x0p = rij.*cos(omega); % (eqn 27)
y0p = rij.*sin(omega); % (eqn 28)

S = sqrt((Xj(2) - Xj(1)).^2 + (Yj(2) - Yj(1)).^2); % plate length
a =-S/2;
b = S/2;

psi=(g_j./(2*pi)).* ( ...
     ((x0p-b)./2).*(log((x0p-b).^2 + y0p.^2)) + y0p.*atan((x0p-b)./y0p)+b ...
   -(((x0p-a)./2).*(log((x0p-a).^2 + y0p.^2)) + y0p.*atan((x0p-a)./y0p)+a));

psi_inf_local =  Xmi*Vinf - Ymi*Uinf;
end


function [ vi,ui,Vinf_i,Uinf_i ] = GET_induced_norm_coeff( Xi, Yi, Xj, Yj, Uinf )
%% DESRIPTION
%  takes the location of two plates i and j and returns the coffecicent
%  for the induced flow normal to the centre of plate i from j

%% IN
%       start and end ploints for two plates
%       Uing ambient flow velocity

%% OUT
%       vi = coefficient for induced flow normal to the centre on plate i
%       from plate j
%
%       ui = coefficient for induced flow tangential to the centre on plate
%       i from plate j
%
%       Vinf_i = component of ambient flow flow normal to the centre of the plate
%
%       Uinf_i = component of ambient flow flow tangential to the centre of the plate 

%set source strength to unity to isolate coefficicent
q_j = 1; 

% panel i
Xmi=0.5*(Xi(2)+Xi(1)); % midpoint
Ymi=0.5*(Yi(2)+Yi(1));
Phi_i=atan2((Yi(2) -Yi(1)),(Xi(2) - Xi(1))); %plate i angle (eqn 24)

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

vj = (q_j./(2*pi)).*(atan(((S./2)-x0p)./y0p)...
    -atan((-(S./2) - x0p)./y0p)); % eqn(30)

uj = (q_j./(2*pi)).*((-log((y0p.^2+((S.^2)./4)- (S.*x0p)+x0p.^2))./2)...
    + (log((y0p.^2 + ((S.^2)./4) + (S.*x0p) + x0p.^2))./2)); % eqn(29)

vi=uj.*sin(Phi_j-Phi_i)+vj.*cos(Phi_j-Phi_i); % eqn(31)
ui=uj.*cos(Phi_j-Phi_i)-vj.*sin(Phi_j-Phi_i); % eqn(32)

Vinf_i = -Uinf*sin(2*pi - Phi_i); % ambient flow normal to plate
Uinf_i = -Uinf*cos(2*pi - Phi_i); % ambient flow tangential to plate


end


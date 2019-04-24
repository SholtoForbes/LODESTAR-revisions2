function [dz,q,xi,vec_angle,T,D,dm,M0,rho0,P0,T0,P_1_tip,T_1_tip,rho_1_tip,M_1_tip,v_1_tip,heating_rate ] = FirstStageDynamics(z,u,t,phase,interp,Vehicle,Atmosphere,auxdata)
% function [dz,q,phi] = rocketDynamics(z,u,t,phase,scattered)
%This function determines the dynamics of the system as time derivatives,
%from input primals. Called by FirstStageDynamics

% global mach
interp = auxdata.interp;

alt = z(1,:);   %Height
v = z(2,:);   %Velocity
m = z(3,:);   %Mass
gamma = z(4,:);
alpha = z(5,:);
zeta = z(6,:);

dalphadt = z(7,:);

phi = z(8,:);

Throttle = z(10,:);

dalphadt2 = u(1,:); % control is second derivative of AoA with time

dThrottledt = u(2,:);

if isnan(gamma)
    gamma = 1.5708;
end

%%%% Compute gravity from inverse-square law:
rEarth = 6.3674447e6;  %(m) radius of earth
mEarth = 5.9721986e24;  %(kg) mass of earth
G = 6.67e-11; %(Nm^2/kg^2) gravitational constant
g = G*mEarth./((alt+rEarth).^2);

g0 = 9.81;

% if h>=0 & strcmp(phase,'postpitch')
% density = interp1(Atmosphere(:,1),Atmosphere(:,4),h);
% P0 = interp1(Atmosphere(:,1),Atmosphere(:,3),h);
% speedOfSound = interp1(Atmosphere(:,1),Atmosphere(:,5),h);
% else
%     density = interp1(Atmosphere(:,1),Atmosphere(:,4),0);
% P0 = interp1(Atmosphere(:,1),Atmosphere(:,3),0);
% speedOfSound = interp1(Atmosphere(:,1),Atmosphere(:,5),0);
% end

if alt>=0 & strcmp(phase,'postpitch')
rho0 = ppval(interp.rho_spline,alt); % Calculate density using atmospheric data
P0 = ppval(interp.P0_spline, alt);
speedOfSound = ppval(interp.c_spline,alt); % Calculate speed of sound using atmospheric data
T0 = ppval(interp.T0_spline, alt); 
P0 = ppval(interp.P0_spline, alt);

else
    rho0 = ppval(interp.rho_spline,0); % Calculate density using atmospheric data
P0 = ppval(interp.P0_spline, 0);
speedOfSound = ppval(interp.c_spline,0); % Calculate speed of sound using atmospheric data
T0 = ppval(interp.T0_spline, 0); 
P0 = ppval(interp.P0_spline, 0);

end

q = 0.5*rho0.*v.^2;
% h
% P_atm
% Merlin 1C engine 
T = Vehicle.T.SL + (101325 - P0)*Vehicle.T.Mod; % Thrust from Falcon 1 users guide. exit area calculated in SCALING.docx
% T
T = T.*Throttle; % Throttle down



Isp = Vehicle.Isp.SL + (101325 - P0)*Vehicle.Isp.Mod; % linear regression of SL and vacuum Isp. From encyclopaedia astronautica, backed up by falcon 1 users guide
% T
% Isp
% g0
dm = -T./Isp./g0;

M0 = v./speedOfSound;

% interpolate coefficients
% Cd = interp.DragGridded(mach,rad2deg(alpha)) ;
% Cl = interp.LiftGridded(mach,rad2deg(alpha)) ;
Cd = interp.DragGridded(M0,rad2deg(alpha)) + interp.Cd_gridded_Visc1(M0,rad2deg(alpha),alt);
Cl = interp.LiftGridded(M0,rad2deg(alpha)) + interp.Cl_gridded_Visc1(M0,rad2deg(alpha),alt);

Cm = (1-(auxdata.Vehicle.mFuel - (m - auxdata.m1FuelDepleted))./auxdata.Vehicle.mFuel).*interp.MomentGriddedFull(M0,rad2deg(alpha)) + (auxdata.Vehicle.mFuel - (m - auxdata.m1FuelDepleted))./auxdata.Vehicle.mFuel.*interp.MomentGriddedEmpty(M0,rad2deg(alpha));

% Compute Thrust Vector

CG = (1-(auxdata.Vehicle.mFuel - (m - auxdata.m1FuelDepleted))./auxdata.Vehicle.mFuel).*auxdata.Vehicle.CG_Full + (auxdata.Vehicle.mFuel - (m - auxdata.m1FuelDepleted))./auxdata.Vehicle.mFuel.*auxdata.Vehicle.CG_Empty;

% vec_angle = -asin(Cm./(T/(CG - auxdata.Vehicle.L-1)./q));

vec_angle = -asin(Cm.*auxdata.Vehicle.L.*auxdata.Vehicle.Area.*q ./(T.*(auxdata.Vehicle.L-CG-1)));
aero_pitching_moment = Cm.*auxdata.Vehicle.L.*auxdata.Vehicle.Area.*q;
thrust_pitching_moment = T.*(auxdata.Vehicle.L-CG-1);
vec_angle(aero_pitching_moment>thrust_pitching_moment) = -deg2rad(89.9);
vec_angle(aero_pitching_moment<-thrust_pitching_moment) = deg2rad(89.9);


%%%% Compute the drag and lift:

Area = Vehicle.Area; 

D = 0.5*Cd.*Area.*rho0.*v.^2*auxdata.Cdmod;
global L
L = 0.5*Cl.*Area.*rho0.*v.^2;

%modify lift and thrust for vectoring
% L = L + T.*sin(vec_angle);
% T = T.*cos(vec_angle);

switch phase
    case 'prepitch'
    gamma = 1.5708*ones(1,length(alt)); % Control Trajectory Angle 
    zeta = 0;
    case 'postpitch'
    %Do nothing
end

xi = zeros(1,length(alt)); 

% This determines the dynamics of the system.
% Set up like this because xi is a quasi-forward sim instead of a primal




% [dr,dxi,dphi,dgamma,dv,dzeta] = RotCoordsFirst(h+rEarth,xi,phi,gamma,v,zeta,L,D,T,m,alpha,phase);
[dr,dxi,dphi,dgamma,dv,dzeta] = RotCoordsFirst(alt,xi,phi,gamma,v,zeta,L,D,T,m,alpha,phase,vec_angle);



switch phase
    case 'prepitch'
    dgamma = 0; % Control Trajectory Angle 
    dzeta = 0;
    case 'postpitch'
    %Do nothing
end

dz = [dr;dv;dm;dgamma;dalphadt;dzeta;dalphadt2;dphi;dxi;dThrottledt];

if any(isnan(dz))
    disp('NaN Values Detected')
end



%% heating---------------------------

% properties after the normal shock at tip
gamma = 1.40; 

P_1_tip(M0>1) = P0(M0>1).*(2*gamma.*M0(M0>1).^2 - (gamma-1))./(gamma+1);

T_1_tip(M0>1) = T0(M0>1).*(2.*gamma.*M0(M0>1).^2 - (gamma-1)).*((gamma-1).*M0(M0>1).^2 + 2)./((gamma+1).^2.*M0(M0>1).^2);

rho_1_tip(M0>1) = rho0(M0>1).*((gamma+1).*M0(M0>1).^2)./((gamma-1).*M0(M0>1).^2+2);

M_1_tip(M0>1) = sqrt(((gamma-1).*M0(M0>1).^2 + 2)./(2*gamma.*M0(M0>1).^2 - (gamma-1)));

P_1_tip(M0<=1) = P0(M0<=1);

T_1_tip(M0<=1) = T0(M0<=1);

rho_1_tip(M0<=1) = rho0(M0<=1);

M_1_tip(M0<=1) = M0(M0<=1);


R = 287.035;
c_1_tip = sqrt(gamma*R*T_1_tip);

v_1_tip = M_1_tip.*c_1_tip;

% Calculate heating rate at tip
% From Conceptual Shape Optimization of Entry Vehicles, Dirkx & Mooj & NASA
% lecture

%using hot wall correction

kappa = 1.83e-4; % sutton-graves, from nasa lecture
Rn = 0.005; %effective nose radius (m) 

% heating_rate = kappa*sqrt(rho0./Rn).*v.^3; %W/m^2

heating_rate = kappa*sqrt(rho_1_tip./Rn).*v_1_tip.^3; %W/m^2



end
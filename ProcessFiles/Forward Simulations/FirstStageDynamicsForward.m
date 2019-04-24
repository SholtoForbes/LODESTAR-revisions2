function [dz] = FirstStageDynamicsForward(z,zeta,alpha,phase,interp,Throttle,Vehicle,Atmosphere,auxdata)
global mach

h = z(1,:);   %Height
v = z(2,:);   %Velocity
m = z(3,:) ;  %Mass
gamma = z(4,:);
phi = z(5,:);


dalphadt = 0;

if isnan(gamma)
    gamma = 1.5708;
end

%%%% Compute gravity from inverse-square law:
rEarth = 6.3674447e6;  %(m) radius of earth
mEarth = 5.9721986e24;  %(kg) mass of earth
G = 6.67e-11; %(Nm^2/kg^2) gravitational constant
g = G*mEarth./((h+rEarth).^2);
if h<80000
density = interp1(Atmosphere(:,1),Atmosphere(:,4),h);
P_atm = interp1(Atmosphere(:,1),Atmosphere(:,3),h);
speedOfSound = interp1(Atmosphere(:,1),Atmosphere(:,5),h);
else
density = 0;
P_atm = 0;
speedOfSound = 0;
end
q = 0.5*density.*v.^2;

SCALE = 1.;
% SCALE = 1; %this is engine exit area scale
% Merlin 1C engine 

T = Vehicle.T.SL + (101325 - P_atm)*Vehicle.T.Mod; % Thrust from Falcon 1 users guide. exit area calculated in SCALING.docx

T = T.*Throttle; % Throttle down

Isp = Vehicle.Isp.SL + (101325 - P_atm)*Vehicle.Isp.Mod; % linear regression of SL and vacuum Isp. From encyclopaedia astronautica, backed up by falcon 1 users guide


dm = -T./Isp./g*SCALE;


mach = v./speedOfSound;

% Cd = interp.DragGridded(mach,rad2deg(alpha));
% Cl = interp.LiftGridded(mach,rad2deg(alpha));

Cd = interp.DragGridded(mach,rad2deg(alpha)) + interp.Cd_gridded_Visc1(mach,rad2deg(alpha),h);
Cl = interp.LiftGridded(mach,rad2deg(alpha)) + interp.Cl_gridded_Visc1(mach,rad2deg(alpha),h);
Cm = (1-(auxdata.Vehicle.mFuel - (m - auxdata.m1FuelDepleted))./auxdata.Vehicle.mFuel).*interp.MomentGriddedFull(mach,rad2deg(alpha)) + (auxdata.Vehicle.mFuel - (m - auxdata.m1FuelDepleted))./auxdata.Vehicle.mFuel.*interp.MomentGriddedEmpty(mach,rad2deg(alpha));

% Compute Thrust Vector
CG = (1-(auxdata.Vehicle.mFuel - (m - auxdata.m1FuelDepleted))./auxdata.Vehicle.mFuel).*auxdata.Vehicle.CG_Full + (auxdata.Vehicle.mFuel - (m - auxdata.m1FuelDepleted))./auxdata.Vehicle.mFuel.*auxdata.Vehicle.CG_Empty;

vec_angle = -asin(Cm.*auxdata.Vehicle.L.*auxdata.Vehicle.Area.*q ./(T.*(auxdata.Vehicle.L-CG-1)));
aero_pitching_moment = Cm.*auxdata.Vehicle.L.*auxdata.Vehicle.Area.*q;
thrust_pitching_moment = T.*(auxdata.Vehicle.L-CG-1);
vec_angle(aero_pitching_moment>thrust_pitching_moment) = -deg2rad(89.9);
vec_angle(aero_pitching_moment<-thrust_pitching_moment) = deg2rad(89.9);
%%%% Compute the drag:

Area = Vehicle.Area; 
D = 0.5*Cd.*Area.*density.*v.^2*auxdata.Cdmod;
global L
L = 0.5*Cl.*Area.*density.*v.^2;

%modify lift and thrust for vectoring
% L = L + T.*sin(vec_angle);
% T = T.*cos(vec_angle);

switch phase
    case 'prepitch'
    gamma = 1.5708*ones(1,length(h)); % Control Trajectory Angle 
    case 'postpitch'
    %Do nothing
end

xi = 0; 

% [dr,dxi,dphi,dgamma,dv,dzeta] = RotCoordsFirst(h+rEarth,xi,phi,gamma,v,zeta,L,D,T,m,alpha,phase);
[dr,dxi,dphi,dgamma,dv,dzeta] = RotCoordsFirst(h,xi,phi,gamma,v,zeta,L,D,T,m,alpha,phase,vec_angle);

% dzeta
if isnan(dgamma)
dgamma = 0;
end

dz = [dr;dv;dm;dgamma;dphi;dzeta];

if any(isnan(dz))
    disp('NaN Values Detected')
end
% dz = [dr;dv;dm;dgamma];
end
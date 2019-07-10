function [rdot,xidot,phidot,gammadot,vdot,zetadot, mdot, Vec_angle, T, L, D, q, Kn] = ThirdStageDynamics(alt,gamma,v,m,Alpha,time,auxdata, Alphadot, phi, zeta)
% Function for simulating the Third Stage Rocket Trajectory
% Created by Sholto Forbes-Spyratos
mode = auxdata.mode;

% mHS = 130.9; % Heat Shield Mass. BASELINE
mHS = 124.6; % Heat Shield Mass

% Calculate approximate fuel mass for CG variation.
% mFuel = m - (auxdata.Stage3.mTot - mHS)*0.09 - mHS;
% mFuel_Full = auxdata.Stage3.mTot - (auxdata.Stage3.mTot - mHS)*0.09 - mHS;


% Set third stage length (m)
L_ThirdStage_ref = 9; % Reference length for cart3d, set to same as clic
L_ThirdStage = 8.66; % Actual Length

time1 = cputime;

% In case of altitude nan inputs, which GPOPS-2 can use, set altitude to 40km
alt(isnan(alt)) = 40000;

Atmosphere = auxdata.interp.Atmosphere;

% Drag_interp = auxdata.interp.Drag_interp3;
% 
% Lift_interp = auxdata.interp.Lift_interp3;
% 
% CP_interp = auxdata.interp.CP_interp3;
% 
% CN_interp = auxdata.interp.CN_interp3;

% Max_AoA_interp = auxdata.interp.Max_AoA_interp3;


% Initiate iterations
iteration = 1;


rho_init = ppval(auxdata.interp.rho_spline, alt(1));
c_init = ppval(auxdata.interp.c_spline, alt(1));

q_init = 0.5.*rho_init.*v.^2;
M_init = v./c_init;

%% Set Max AoA
% determine the maximum allowable normal coefficient with a 10 degree limit
% at 50kPa dynamic pressure, and set the max AoA to match this normal
% coefficient
Alt_50 = spline( Atmosphere(:,4),  Atmosphere(:,1), 50000.*2./v(1).^2);
c_50 = ppval(auxdata.interp.c_spline, Alt_50);

M_50 = v(1)./c_50;
% M_50 = 2922.8./c_50;% using a constant velocity

% CN_50 = CN_interp(M_50,10, Alt_50); % note this is not modified for hypercube, because both the cn at 50, and the max aoa cn will be modified and cancel out


% CN_50 = CN_interp(M_50,5);
% AoA_max = deg2rad(Max_AoA_interp(M_init,CN_50.*50000./q_init)); %maximum allowable AoA

%Reference area

%BASELINE
A = 0.95; % diameter of 1.1m BASELINE

% % Diameter 0.9m
% A = 0.63617;

% Diameter 1m
% A = 0.785;

g = 9.806; %standard gravity

% the Isp influences the optimal burn mass
% Isp = 437; % from Tom Furgusens Thesis %RL10

% BASELINE
if mode ~= 1000
Isp = 317.*0.98* auxdata.Isp3mod; %Kestrel, from Falcon 1 users guide, with efficiency reduction from area reduction. BASELINE
else
Isp = 317.*0.98 ; % if mode = 1000 isp3mod is used later
end
% % Diameter 0.9m
% Isp = 317.*0.96; %Kestrel, from Falcon 1 users guide, with efficiency reduction.

% Diameter 1m
% Isp = 317.*0.97; %Kestrel, from Falcon 1 users guide, with efficiency reduction. 

% Isp = 317; % Kestrel
% Isp = 446; %HM7B
% Isp = 340; %Aestus 2

%% Define Vehicle Properties
%BASELINE

% %Diameter 0.9m
% mHS = 109.3;

%Diameter 1m
% mHS = 120.341;

% mEng = 100; %RL10
% mEng = 52; %Kestrel. BASELINE
% mEng = 165; %HM7B
% mEng = 138; %Aestus 2 ./ RS72 from https:././web.archive.org./web./20141122143945./http:././cs.astrium.eads.net./sp./launcher-propulsion./rocket-engines./aestus-rs72-rocket-engine.html
mEng = 78; %Kestrel. Mass Increased.

% mdot = 14.71; %RL10
% mdot = 9.86977; %Kestrel
mdot = 9.86977.*1.5; %Kestrel Modified
% mdot = 9.86977.*1.4; %Kestrel Modified
% mdot = 14.8105; %HM7B
% mdot = 16.5; %Aestus 2

%BASELINE
% Calculate approximate fuel mass for CG variation.
mFuel = m - (auxdata.Stage3.mTot - mHS - mEng - 100)*0.079 - mHS;
mFuel_Full = auxdata.Stage3.mTot - (auxdata.Stage3.mTot - mHS - mEng - 100)*0.079 - mHS;

% CG
CG_Full = 4.326; %CG Full of fuel. 

CG_Empty = 4.046; % CG no fuel

CG = mFuel/mFuel_Full.*CG_Full + (1-mFuel/mFuel_Full).*CG_Empty;


% % Diameter 0.9m
% CG =  5.110-0.7521;

% Set the moment of inertia
%Baseline
I = 6590; % From Creo

% % Diameter 0.9m
% I = 4864;

% Diameter 1m
% I = 5700;

%% Calculate Aerodynamics

% Define the altitude that the rocket exits the atmosphere
atmo_elem = find(alt<=85000); % Find the elements that are within the atmosphere
exo_elem = find(alt>85000); % Find the elements that are exoatmospheric

% Calculate exoatmospheric properties
rho(atmo_elem) = ppval(auxdata.interp.rho_spline, alt(atmo_elem)); % Density kg/m^3, in atmosphere
c(atmo_elem) = ppval(auxdata.interp.c_spline, alt(atmo_elem)); % Speed of sound (m/s), in atmosphere
p(atmo_elem) = ppval(auxdata.interp.P0_spline, alt(atmo_elem)); % Pressure (Pa), in atmosphere
T(atmo_elem) = ppval(auxdata.interp.T0_spline, alt(atmo_elem)); % Temperatue (K), in atmosphere

% Calculate exoatmospheric properties, outside of atmosphere
% This smooths the atmospheric properties to 0 as the rocket leaves the
% atmosphere
rho(exo_elem) = ppval(auxdata.interp.rho_spline,85000).*gaussmf(alt(exo_elem),[100 85000]); % Density kg/m^3, exoatmosphere
c(exo_elem) = ppval(auxdata.interp.c_spline, 85000); % Speed of sound (m/s), exoatmosphere
p(exo_elem) = ppval(auxdata.interp.p_spline,85000).*gaussmf(alt(exo_elem),[100 85000]); % Pressure (Pa), exoatmosphere
T(exo_elem) = ppval(auxdata.interp.T0_spline, alt(exo_elem)); % Temperatue (K), in atmosphere


rho = rho.';
c = c.';
p = p.';

q(atmo_elem) = 0.5.*rho(atmo_elem).*v(atmo_elem).^2; % Dynamic pressure (Pa), in atmosphere
M(atmo_elem)  = v(atmo_elem) ./c(atmo_elem) ; % Mach number, in atmosphere
% AoA_max(atmo_elem) = deg2rad(Max_AoA_interp(M(atmo_elem),CN_50.*50000./q(atmo_elem))); %maximum allowable AoA, in atmosphere

q(exo_elem) = 0.5.*rho(exo_elem).*v(exo_elem).^2; % Dynamic pressure (Pa), exoatmosphere
M(exo_elem)  = v(exo_elem) ./c(exo_elem) ; % Mach number, exoatmosphere
% AoA_max(exo_elem) = deg2rad(30); %maximum allowable AoA, exoatmosphere

q = q.';
M = M.';
% AoA_max = AoA_max.';

% Calculate thrust, assuming that the exit area is equal to the total area
% of the rocket

T  = Isp.*mdot.*g - p .*A; % Thrust (N)
% T = T* auxdata.Isp3mod; % Modify Thrust if necessary

% Interpolate for aerodynamic coefficients
% CD  = Drag_interp(M ,rad2deg(Alpha ),alt);
% CL  = Lift_interp(M ,rad2deg(Alpha ),alt);
% CN  = CN_interp(M ,rad2deg(Alpha ),alt);
% cP  = CP_interp(M ,rad2deg(Alpha ),alt);




CD  = auxdata.interp.D3interp(M ,rad2deg(Alpha ),alt);
CL  = auxdata.interp.L3interp(M ,rad2deg(Alpha ),alt);


if mode == 1000

CD  = CD + CD*auxdata.Cn3mod.*sin(Alpha)+ CD*auxdata.Ca3mod.*cos(Alpha);
CL  = CL + CL*auxdata.Cn3mod.*cos(Alpha)+ CL*auxdata.Ca3mod.*sin(Alpha);
CN  = CN + CN*auxdata.Cm3mod; % This is modified by the moment modifier because it is only used to calculate moment, so this modification serves to modify the moment directly
T = T + T*auxdata.Isp3mod;

end

% Calculate aerodynamic forces
D  = 1./2.*rho.*(v.^2).*A.*CD *auxdata.Cd3mod;
L  = 1./2.*rho.*(v.^2).*A.*CL ; % Aerodynamic lift
% N  = 1./2.*rho.*(v.^2).*A.*CN ;
    
   
%% Thrust vectoring
% 1.5m is subtracted as the CG is taken from the end of the nozzle, to make
% the cg relative to the base of the rocket fuselage, where the thrust
% force is applied

%cP is from nose, negative along body length

% Vec_angle  = asin((CG+(cP.*1.1))./(L_ThirdStage-CG-1.5).*(N)./T - Alphadot.*I./((L_ThirdStage-CG-1.5).*T)); % calculate the thrust vector angle necessary to resist the lift force moment. cP is in ref lengths. note cP is negative from tip

Cm_full = auxdata.interp.momentInterp3_full(M ,rad2deg(Alpha ));
Cm_empty = auxdata.interp.momentInterp3_empty(M ,rad2deg(Alpha ));

Cm = (CG-CG_Empty)./(CG_Full-CG_Empty).*Cm_full + (CG_Full - CG)./(CG_Full-CG_Empty).*Cm_empty;


Vec_angle  = -asin(Cm.*L_ThirdStage_ref.*A.*q ./(T.*(L_ThirdStage_ref-CG-1.5)));

if any(not(isreal(Vec_angle)))
Vec_angle(not(isreal(Vec_angle)))  = deg2rad(80); % This stops the vector angle going imaginary
end


%% Calculate Dynamics
% [rdot,xidot,phidot,gammadot,vdot,zetadot] = RotCoordsRocket(alt+auxdata.Re,0,phi,gamma,v,zeta,L,D,T,m,Alpha,Vec_angle);
% %Also change this in third stage sim
[rdot,xidot,phidot,gammadot,vdot,zetadot] = RotCoordsRocket(alt,zeros(length(alt),1),phi,gamma,v,zeta,L,D,T,m,Alpha,Vec_angle);

end


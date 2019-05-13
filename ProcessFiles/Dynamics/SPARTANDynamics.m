function [altdot,xidot,phidot,gammadot,a,zetadot, q, M0, D, rho0,L,Fueldt,T,Isp,q1,flap_deflection,heating_rate_stag,CG,T1,P1,M1,P0,T0,P_1_tip,T_1_tip,rho_1_tip,M_1_tip,v_1_tip,heating_rate_LE] = SPARTANDynamics(gamma, alt, v,auxdata,zeta,phi,xi,alpha,eta,throttle,mFuel,mFuelinit,mFuelend,ThirdStage,forwardflag)
%===================================================
%
% SPARTAN DYNAMICS SIMULATION
%
%===================================================

interp = auxdata.interp;
mode = auxdata.mode;
% =======================================================
% Vehicle Model
% =======================================================
A = auxdata.A; % reference area (m^2)

if ThirdStage == 1
m = auxdata.Stage2.mStruct+mFuel+auxdata.Stage3.mTot; 
else
m = auxdata.Stage2.mStruct+mFuel;
end

%% Flow =============================================================
c = ppval(interp.c_spline,alt); % Calculate speed of sound using atmospheric data
mach = v./c;
rho0 = ppval(interp.rho_spline,alt); % Calculate density using atmospheric data

q = 0.5 * rho0 .* (v .^2); % Calculating Dynamic Pressure

M0 = v./c; % Calculating Mach No (Descaled)

T0 = ppval(interp.T0_spline, alt); 

P0 = ppval(interp.P0_spline, alt);
%% Thrust 


[Isp_nozzlefront,Fueldt_max,eq,q1,T1,P1,M1] = RESTint(M0, rad2deg(alpha), auxdata,T0,P0); % Calculate C-REST engine properties

% Isp = Isp_nozzlefront.*auxdata.Ispmod ; %
Isp = Isp_nozzlefront; %

% Isp(q1<20000) = Isp(q1<20000).*gaussmf(q1(q1<20000),[1000,20000]); % rapidly reduce ISP to 0 after passing the lower limit of 20kPa dynamic pressure. This dynamic pressure is after the conical shock.


 if ThirdStage == 0 && forwardflag ==0

    % Turn off throttle at unoperable flight conditions for aerodynamic and
    % thrust purposes
    throttle(q1<20000) = throttle(q1<20000).*gaussmf(q1(q1<20000),[100,20000]); % rapidly reduce throttle to 0 after passing the lower limit of 20kPa dynamic pressure. This dynamic pressure is after the conical shock.
    throttle(M0<5.0) =   throttle(M0<5.0).*gaussmf(M0(M0<5.0),[.01,5]); % remove throttle points below operable range on return flight
 else
     throttle(q1<20000) = throttle(q1<20000).*gaussmf(q1(q1<20000),[100,20000]);
end

Fueldt = Fueldt_max.*throttle; %

% T = Isp.*Fueldt*9.81.*cos(alpha).*gaussmf(throttle,[0.1,1]); % Thrust in
% direction of , modified by a gaussmf funtion to reduce thrust rapidly

if auxdata.mode == 3
T = Isp.*Fueldt_max.*throttle*9.81 + auxdata.T_spline(mach,rad2deg(alpha),alt/1000)*(auxdata.Ispmod-1).*throttle; % Thrust in direction of motion, if Isp is modified, add portion of total thrust
elseif mode ~= 1000
T = Isp.*Fueldt_max.*throttle*9.81;
end
%% Aerodynamics
% Calculate aerodynamic coefficients

   
    % Interpolate between centre of gravity conditions for full, cylindrical tank empty, and empty conditions as fuel depletes
mFuel_cyltanks = 710;

CG_withFuel_noThirdStage = 14.55; % CG from CREO with no third stage but full fuel
CG_cyltanksEmpty_noThirdStage = 14.30;
CG_cyltanksEmpty_withThirdStage = 15.13; % CG with third stage, but no fuel
CG_noFuel_noThirdStage = 15.16;% CG  at no fuel conditition (it is assumed that the fuel for the return is used so as to not change the CG)
CG_noFuel_withThirdStage = 15.74; % CG with third stage, but no fuel
CG_withFuel_withThirdStage = 15.24; % CG with third stage and full fuel
% 

if ThirdStage == 1 % Ascent with third stage

% Determine trajectory points which are using fuel from cylindrical fuel tank
index_overcyl = mFuel>(auxdata.Stage2.mFuel-mFuel_cyltanks);

% Determine the fraction of fuel in the cylindrical tank
Proportion_fulltocyl  = (mFuel(index_overcyl)-(auxdata.Stage2.mFuel-mFuel_cyltanks))./(mFuel_cyltanks);

% Calculate CG variation
CG(index_overcyl) = Proportion_fulltocyl*CG_withFuel_withThirdStage + (1-Proportion_fulltocyl)*CG_cyltanksEmpty_withThirdStage;


% Calculate Cd
Cd(index_overcyl) = Proportion_fulltocyl.*auxdata.interp.Cd_spline_EngineOn.fullFuel(mach(index_overcyl),...
    rad2deg(alpha(index_overcyl)),alt(index_overcyl)/1000) + (1-Proportion_fulltocyl).*...
    auxdata.interp.Cd_spline_EngineOn.cylTankEnd(mach(index_overcyl),rad2deg(alpha(index_overcyl)),alt(index_overcyl)/1000);

% Calculate Cl
Cl(index_overcyl) = Proportion_fulltocyl.*auxdata.interp.Cl_spline_EngineOn.fullFuel(mach(index_overcyl),...
    rad2deg(alpha(index_overcyl)),alt(index_overcyl)/1000) + (1-Proportion_fulltocyl).*...
    auxdata.interp.Cl_spline_EngineOn.cylTankEnd(mach(index_overcyl),rad2deg(alpha(index_overcyl)),alt(index_overcyl)/1000);

% Calculate frap deflection
flap_deflection(index_overcyl) = Proportion_fulltocyl.*auxdata.interp.flap_spline_EngineOn.fullFuel(mach(index_overcyl),...
    rad2deg(alpha(index_overcyl)),alt(index_overcyl)/1000) + (1-Proportion_fulltocyl).*...
    auxdata.interp.flap_spline_EngineOn.cylTankEnd(mach(index_overcyl),rad2deg(alpha(index_overcyl)),alt(index_overcyl)/1000);





% Determine trajectory points which are using fuel from the conical tank
index_undercyl = mFuel<=(auxdata.Stage2.mFuel-mFuel_cyltanks);

% Determine fraction of fuel in conical tank
Proportion_EmptytoCyl = ((auxdata.Stage2.mFuel-mFuel_cyltanks)-mFuel(index_undercyl))./(auxdata.Stage2.mFuel-mFuel_cyltanks);

% Calculate CG
CG(index_undercyl) = Proportion_EmptytoCyl*CG_noFuel_withThirdStage + (1-Proportion_EmptytoCyl)*CG_cyltanksEmpty_withThirdStage;

Cd(index_undercyl) = Proportion_EmptytoCyl.*auxdata.interp.Cd_spline_EngineOn.noFuel(mach(index_undercyl),...
    rad2deg(alpha(index_undercyl)),alt(index_undercyl)/1000) + (1-Proportion_EmptytoCyl).*...
    auxdata.interp.Cd_spline_EngineOn.cylTankEnd(mach(index_undercyl),rad2deg(alpha(index_undercyl)),alt(index_undercyl)/1000);

Cl(index_undercyl) = Proportion_EmptytoCyl.*auxdata.interp.Cl_spline_EngineOn.noFuel(mach(index_undercyl),...
    rad2deg(alpha(index_undercyl)),alt(index_undercyl)/1000) + (1-Proportion_EmptytoCyl).*...
    auxdata.interp.Cl_spline_EngineOn.cylTankEnd(mach(index_undercyl),rad2deg(alpha(index_undercyl)),alt(index_undercyl)/1000);

flap_deflection(index_undercyl) = Proportion_EmptytoCyl.*auxdata.interp.flap_spline_EngineOn.noFuel(mach(index_undercyl),...
    rad2deg(alpha(index_undercyl)),alt(index_undercyl)/1000) + (1-Proportion_EmptytoCyl).*...
auxdata.interp.flap_spline_EngineOn.cylTankEnd(mach(index_undercyl),rad2deg(alpha(index_undercyl)),alt(index_undercyl)/1000);

Cd = Cd.';

Cl = Cl.' ;
flap_deflection = flap_deflection.';

else  % Return, no third stage
    
    
Proportion_EmptytoCyl = ((auxdata.Stage2.mFuel-mFuel_cyltanks)-mFuel)./(auxdata.Stage2.mFuel-mFuel_cyltanks);

CG = Proportion_EmptytoCyl*CG_noFuel_noThirdStage + (1-Proportion_EmptytoCyl)*CG_cyltanksEmpty_noThirdStage;


Cd_Engineoff = Proportion_EmptytoCyl.*auxdata.interp.Cd_spline_EngineOff.noFuel(mach,...
    rad2deg(alpha),alt/1000) + (1-Proportion_EmptytoCyl).*...
    auxdata.interp.Cd_spline_EngineOff.cylTankEnd(mach,rad2deg(alpha),alt/1000);

Cl_Engineoff = Proportion_EmptytoCyl.*auxdata.interp.Cl_spline_EngineOff.noFuel(mach,...
    rad2deg(alpha),alt/1000) + (1-Proportion_EmptytoCyl).*...
    auxdata.interp.Cl_spline_EngineOff.cylTankEnd(mach,rad2deg(alpha),alt/1000);

flap_deflection_Engineoff = Proportion_EmptytoCyl.*auxdata.interp.flap_spline_EngineOff.noFuel(mach,...
    rad2deg(alpha),alt/1000) + (1-Proportion_EmptytoCyl).*...
auxdata.interp.flap_spline_EngineOff.cylTankEnd(mach,rad2deg(alpha),alt/1000);

Cd_Engineon = Proportion_EmptytoCyl.*auxdata.interp.Cd_spline_EngineOn.noFuel(mach,...
    rad2deg(alpha),alt/1000) + (1-Proportion_EmptytoCyl).*...
    auxdata.interp.Cd_spline_EngineOn.cylTankEnd(mach,rad2deg(alpha),alt/1000);

Cl_Engineon = Proportion_EmptytoCyl.*auxdata.interp.Cl_spline_EngineOn.noFuel(mach,...
    rad2deg(alpha),alt/1000) + (1-Proportion_EmptytoCyl).*...
    auxdata.interp.Cl_spline_EngineOn.cylTankEnd(mach,rad2deg(alpha),alt/1000);

flap_deflection_Engineon = Proportion_EmptytoCyl.*auxdata.interp.flap_spline_EngineOn.noFuel(mach,...
    rad2deg(alpha),alt/1000) + (1-Proportion_EmptytoCyl).*...
auxdata.interp.flap_spline_EngineOn.cylTankEnd(mach,rad2deg(alpha),alt/1000);

 Cd = (1-throttle).*Cd_Engineoff + throttle.*Cd_Engineon;
Cl = (1-throttle).*Cl_Engineoff + throttle.*Cl_Engineon;
flap_deflection = (1-throttle).*flap_deflection_Engineoff + throttle.*flap_deflection_Engineon;
   
    
end

if ThirdStage == 1 && auxdata.mode == 5
    
Cd = Cd + auxdata.Cd_spline_ViscousEngineOn(mach,rad2deg(alpha),alt/1000).*(auxdata.vCdmod-1);

elseif ThirdStage == 0 && auxdata.mode == 5
        
vCd_add = (1-throttle).*auxdata.Cd_spline_ViscousEngineOff(mach,rad2deg(alpha),alt/1000).*(auxdata.vCdmod-1)...
    + throttle.*auxdata.Cd_spline_ViscousEngineOn(mach,rad2deg(alpha),alt/1000).*(auxdata.vCdmod-1);  
Cd = Cd + vCd_add;
end

if mode == 1000
    % Modify coefficients by modifiers, using sigma mf to smoothly
    % transition between regimes. Regime switches at M=1, with a smooth
    % transition

    Cl(M0<1) = Cl(M0<1) + Cl(M0<1).*(auxdata.CL12_subsonicmod.*(1-sigmf(M0(M0<1),[100,0.8])) + auxdata.CL12_transonicmod.*sigmf(M0(M0<1),[100,0.8]) + auxdata.CL12_supersonicmod.*sigmf(M0(M0<1),[100,1.2]));
    Cl(M0>=1) = Cl(M0>=1) + Cl(M0>=1).*(auxdata.CL12_subsonicmod.*(1-sigmf(M0(M0>=1),[100,0.8])) + auxdata.CL12_transonicmod.*(1-sigmf(M0(M0>=1),[100,1.2])) + auxdata.CL12_supersonicmod.*sigmf(M0(M0>=1),[100,1.2]));
    
    Cd(M0<1) = Cd(M0<1) + Cd(M0<1).*(auxdata.CD12_subsonicmod.*(1-sigmf(M0(M0<1),[100,0.8])) + auxdata.CD12_transonicmod.*sigmf(M0(M0<1),[100,0.8]) + auxdata.CD12_supersonicmod.*sigmf(M0(M0<1),[100,1.2]));
    Cd(M0>=1) = Cd(M0>=1) + Cd(M0>=1).*(auxdata.CD12_subsonicmod.*(1-sigmf(M0(M0>=1),[100,0.8])) + auxdata.CD12_transonicmod.*(1-sigmf(M0(M0>=1),[100,1.2])) + auxdata.CD12_supersonicmod.*sigmf(M0(M0>=1),[100,1.2]));

    T = Isp.*Fueldt_max.*throttle*9.81 + auxdata.T_spline(mach,rad2deg(alpha),alt/1000)*(auxdata.Isp2mod).*throttle; % If Isp is modified, add portion of total thrust
end

%%%% Compute the drag and lift:
D = 0.5*Cd.*A.*rho0.*v.^2*auxdata.Cdmod;
L = 0.5*Cl.*A.*rho0.*v.^2;

%Rotational Coordinates =================================================



[altdot,xidot,phidot,gammadot,a,zetadot] = RotCoords(alt,xi,phi,gamma,v,zeta,L,D,T,m,alpha,eta,auxdata.delta);

% Aero Data =============================================================

q = 0.5 * rho0 .* (v .^2); % Calculating Dynamic Pressure

M0 = v./c; % Calculating Mach No (Descaled)

v_H = v.*cos(gamma);

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

kappa = 1.83e-4; % sutton-graves, from nDirkx
Rn = 0.05; %effective nose radius (m) 

% heating_rate = kappa*sqrt(rho0./Rn).*v.^3; %W/m^2

% heating_rate_stag = kappa*sqrt(rho_1_tip./Rn).*v_1_tip.^3; %W/m^2
heating_rate_stag = kappa*sqrt(rho0./Rn).*v.^3; %W/m^2
% heating_rate_stag = heating_rate_stag';


% Heat Transfer at Wing Leading Edge
sweep_angle = 72.9; % deg

% Laminar
x_l = 5.72; %m, half way along wing, running variable, not sure if this is correct
k_FP = 2.53e-5*cos(deg2rad(sweep_angle)).^0.5*sin(deg2rad(sweep_angle))*x_l.^-0.5; % Assumes cold wall

heating_rate_FP = k_FP*(rho_1_tip./Rn).^0.5.*v_1_tip.^3.2; % laminar powers

heating_rate_FP = heating_rate_FP';

R_wing = 0.0015; % tip radius of wing
heating_rate_stag_wing = kappa*sqrt(rho_1_tip./R_wing).*v_1_tip.^3; %W/m^2
heating_rate_stag_wing = heating_rate_stag_wing';

heating_rate_LE = (0.5*heating_rate_stag_wing*cos(deg2rad(sweep_angle)).^2 + heating_rate_FP*sin(deg2rad(sweep_angle)).^2) ;

% =========================================================================
end









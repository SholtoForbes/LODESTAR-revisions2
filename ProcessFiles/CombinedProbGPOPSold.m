%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rocket-Scramjet-Rocket Launch Optimiser
% By Sholto Forbes-Spyratos
% Utilises the GPOPS-2 proprietary optimisation software
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc

%%
% =========================================================================
% SET RUN MODE
% =========================================================================
% Change mode to set the target of the simulation. Much of the problem
% definition changes with mode.

% mode 1: 50kPa, standard trajectory, used for optimal trajectory calculation
% mode 2: Dynamic Pressure Variation
% mode 3: Isp Variation
% mode 4: Cd Variation
% mode 5: Viscous Cd Variation - goes from approximately fully laminar to
% fully turbulent

% mode 6: First Stage mass
% mode 7: Third Stage Drag Variation
% mode 8: Third Stage Mass Variation
% mode 9: Third Stage Isp variation
% mode 10: SPARTAN mass 
% mode 11: SPARTAN fuel mass 

% mode 12: SPARTAN Thrust angle - not implemented
% mode 13: Boat tail lift - not implemente
% mode 14: Boat tail drag - not implemented

% mode 90: const q
% mode 99: interaction mode

% mode 0: Alternate launch location

mode = 1
auxdata.mode = mode;

returnMode = 1% Flag for setting the return of the SPARTAN. 0 = not constrained (no return), 1 = constrained (return)
auxdata.returnMode = returnMode;




%%

% Set modifiers for mode variation
auxdata.Ispmod = 1;
auxdata.Cdmod = 1;
auxdata.vCdmod = 1;
auxdata.m3mod = 1;
auxdata.Isp3mod = 1;
auxdata.Cd3mod = 1;

addpath('.\Processing\num2words')

% Mode 1
if mode == 1

namelist{1} = 'Standard';
end

% mode 0
if mode == 0

namelist{1} = 'Alternate';
end

if mode == 90

namelist{1} = 'Constq';
end
% Mode 2
q_vars = [40000 45000 50000 55000 60000]  % Set dynamic pressures to be investigated
if mode == 2

%     for i = 1:length(q_vars)
%         namelist{i} = ['qmax' num2words(q_vars(i)/1000)];
%     end
if returnMode == 1
namelist = {'qForty' 'qFortyFive' 'qStandard' 'qFiftyFive' 'qSixty'}
else
namelist = {'qFortyNoReturn' 'qFortyFiveNoReturn' 'qStandardNoReturn' 'qFiftyFiveNoReturn' 'qSixtyNoReturn'}
end
end

% Mode 3
Isp_vars = [0.9 0.95 1 1.05 1.1] 
if mode == 3

%     for i = 1:length(Isp_vars)
%         namelist{i} = ['Isp' num2words(Isp_vars(i)*100) '%'];
%     end
if returnMode == 1
    namelist = {'IspNinety' 'IspNinetyFive' 'IspStandard' 'IspOneHundredFive' 'IspOneHundredTen'}
    else
namelist = {'IspNinetyNoReturn' 'IspNinetyFiveNoReturn' 'IspStandardNoReturn' 'IspOneHundredFiveNoReturn' 'IspOneHundredTenNoReturn'}
end
end
% Mode 4
Cd_vars = [0.9 0.95 1 1.05 1.1] 
if mode == 4
 
%     for i = 1:length(Cd_vars)
%         namelist{i} = ['Cd' num2words(Cd_vars(i)*100) '%'];
%     end 
if returnMode == 1
namelist = {'CdNinety' 'CdNinetyFive' 'CdStandard' 'CdOneHundredFive' 'CdOneHundredTen'}
else
namelist = {'CdNinetyNoReturn' 'CdNinetyFiveNoReturn' 'CdStandardNoReturn' 'CdOneHundredFiveNoReturn' 'CdOneHundredTenNoReturn'}
end
end

% Mode 5
vCd_vars = [0.2 0.5 1 1.07 1.15]
if mode == 5 
%     for i = 1:length(vCd_vars)
%         namelist{i} = ['vCd' num2words(vCd_vars(i)*100) '%'];
%     end  
if returnMode == 1
namelist = {'vCdTwenty' 'vCdFifty' 'vCdStandard' 'vCdOneHundredSeven' 'vCdOneHundredFifteen'}
else
namelist = {'vCdTwentyNoReturn' 'vCdFiftyNoReturn' 'vCdStandardNoReturn' 'vCdOneHundredSevenNoReturn' 'vCdOneHundredFifteenNoReturn'}
end
end
    
% Mode 6
% First Stage Structural Mass (total mass kept constant)
FirstStagem_vars = [0.90 0.95 1 1.05 1.1] 
if mode == 6 
if returnMode == 1
namelist = {'FirstStagemNinety' 'FirstStagemNinetyFive' 'FirstStagemStandard' 'FirstStagemOneHundredFive' 'FirstStagemOneHundredTen'}
else
namelist = {'FirstStagemNinetyNoReturn' 'FirstStagemNinetyFiveNoReturn' 'FirstStagemStandardNoReturn' 'FirstStagemOneHundredFiveNoReturn' 'FirstStagemOneHundredTenNoReturn'}
end
end

% Mode 7
Cd3_vars = [0.80 0.9 1 1.1 1.2]  
if mode == 7
 if returnMode == 1
        namelist = {'CdThreeEighty' 'CdThreeNinety' 'CdThreeStandard' 'CdThreeOneHundredTen' 'CdThreeOneHundredTwenty'};
 else
     namelist = {'CdThreeEightyNoReturn' 'CdThreeNinetyNoReturn' 'CdThreeStandardNoReturn' 'CdThreeOneHundredTenNoReturn' 'CdThreeOneHundredTwentyNoReturn'};
 end
end   
   
% Mode 8
m3_vars = [0.9 0.95 1 1.05 1.1] 
if mode == 8 
%     for i = 1:length(m3_vars)
%         namelist{i} = ['m3' num2words(m3_vars(i)*100) '%'];
%     end  
if returnMode == 1
namelist = {'mThreeNinety' 'mThreeNinetyFive' 'mThreeStandard' 'mThreeOneHundredFive' 'mThreeOneHundredTen'}
else
namelist = {'mThreeNinetyNoReturn' 'mThreeNinetyFiveNoReturn' 'mThreeStandardNoReturn' 'mThreeOneHundredFiveNoReturn' 'mThreeOneHundredTenNoReturn'}
end
end 
    
% Mode 9
Isp3_vars = [0.95 0.975 1 1.025 1.05]  
if mode == 9
%     for i = 1:length(Isp3_vars)
%         namelist{i} = ['T3' num2words(Isp3_vars(i)*100) '%'];
%     end 
if returnMode == 1
namelist = {'ISPThreeNinetyFive' 'ISPThreeNinetySevenFive' 'ISPThreeStandard' 'ISPThreeOneHundredTwoFive' 'ISPThreeOneHundredFive'}
else
namelist = {'ISPThreeNinetyFiveNoReturn' 'ISPThreeNinetySevenFiveNoReturn' 'ISPThreeStandardNoReturn' 'ISPThreeOneHundredTwoFiveNoReturn' 'ISPThreeOneHundredFiveNoReturn'}
end
end  
    
% mode 10
mSPARTAN_vars = [0.95 0.975 1 1.025 1.05] 
if mode == 10 
%     for i = 1:length(mSPARTAN_vars)
%         namelist{i} = ['mSPARTAN' num2words(mSPARTAN_vars(i)*100) '%'];
%     end   
if returnMode == 1
namelist = {'mSPARTANNinetyFive' 'mSPARTANNinetySevenFive' 'mSPARTANStandard' 'mSPARTANOneHundredTwoFive' 'mSPARTANOneHundredFive'}
else
namelist = {'mSPARTANNinetyFiveNoReturn' 'mSPARTANNinetySevenFiveNoReturn' 'mSPARTANStandardNoReturn' 'mSPARTANOneHundredTwoFiveNoReturn' 'mSPARTANOneHundredFiveNoReturn'}
end
end 
 
%mode 11
mFuel_vars = [0.9 0.95 1 1.05 1.1]     
if mode == 11
%     for i = 1:length(mFuel_vars)
%         namelist{i} = ['mFuel' num2words(mFuel_vars(i)*100) '%'];
%     end  
if returnMode == 1
namelist = {'mFuelNinety' 'mFuelNinetyFive' 'mFuelStandard' 'mFuelOneHundredFive' 'mFuelOneHundredTen'}
else
namelist = {'mFuelNinetyNoReturn' 'mFuelNinetyFiveNoReturn' 'mFuelStandardNoReturn' 'mFuelOneHundredFiveNoReturn' 'mFuelOneHundredTenNoReturn'}
end

end

%interactionmode
if mode == 99
 namelist = {};
end

auxdata.namelist = namelist;

%% Misc Modifiers
auxdata.delta = deg2rad(0); % thrust vector angle 


%% Launch Point
% lat0 = deg2rad(-12.4466); % Equatorial Launch Australia Spaceport near Nhulunbuy
% lon0 = deg2rad(136.845);

if mode == 0
% lat0 = deg2rad(-38.4314); % Southern Launch Australia, assumed launch from Cape Nelson
% lon0 = deg2rad(141.5420);
lat0 = deg2rad(-32.9106); % Southern Launch Australia, assumed launch from Streaky Bay
lon0 = deg2rad(134.0724);
%streaky bay council report with support for launch site: https://www.streakybay.sa.gov.au/webdata/resources/minutesAgendas/Council%20Agenda-Report%20-%20October%202017.pdf
else
lat0 = deg2rad(-12.4466); % Equatorial Launch Australia Spaceport near Nhulunbuy
lon0 = deg2rad(136.845);
end


%% Add Necessary Paths


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('..\EngineData')
addpath('..\')
addpath('.\Interpolators')
addpath('.\Processing')
addpath('.\Forward Simulations')
addpath('.\Dynamics')
addpath ('..\ViscousAero')
addpath('..\CG15.1255') % Add folder containing aerodynamic datasets


%% FIRST STAGE


%  auxdata.Throttle = .7 % throttle the Merlin engine down by a modeant value, to enable easier pitchover

% Aerodynamics File Path
Aero = dlmread('FirstStageAero23.5262');


%% Vehicle 

% mRocket =21816; % total mass of scaled Falcon at 9.5m, note, this will not be the final total mass. Calculated using the method outlined in SIZING.docx
mRocket =19569; % total mass of scaled Falcon at 8.5m, note, this will not be the final total mass. Calculated using the method outlined in SIZING.docx

mEngine = 470; % Mass of Merlin 1C
mFuel = 0.939*(mRocket-mEngine); % structural mass fraction calculated without engine
mSpartan = 9819.11;

% Thrust and Isp are modified with altitude through the formula: %
% SL + (101325-P_atm)*Mod %
 
auxdata.Vehicle.T.SL = 555900; % Thrust from Falcon 1 users guide. 
auxdata.Vehicle.T.Mod = 0.5518; % exit area calculated in SIZING.docx

auxdata.Vehicle.Isp.SL = 275; % linear regression of SL and vacuum Isp. From encyclopaedia astronautica, backed up by falcon 1 users guide
auxdata.Vehicle.Isp.Mod = 2.9410e-04;

auxdata.Vehicle.Area = 62.77; 

auxdata.Vehicle.CG = 23.5262; % centre of gravity (m) calculated in clicalcCGvar.m
auxdata.Vehicle.L = 22.94+9.5;
%% Import Atmosphere

auxdata.Atmosphere = dlmread('atmosphere.txt');

%% Calculate Aerodynamic Splines

interp.Lift = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,3));
interp.Drag = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,4));
interp.Moment = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,4));

M_list = unique(sort(Aero(:,1))); % create unique list of Mach numbers from engine data
M_interp = unique(sort(Aero(:,1)));

AoA_list = unique(sort(Aero(:,2))); % create unique list of angle of attack numbers from engine data 
AoA_interp = unique(sort(Aero(:,2)));

[grid.M,grid.AoA] =  ndgrid(M_interp,AoA_interp);
grid.Lift = interp.Lift(grid.M,grid.AoA);
grid.Drag = interp.Drag(grid.M,grid.AoA);
grid.Moment = interp.Moment(grid.M,grid.AoA);

auxdata.interp.LiftGridded = griddedInterpolant(grid.M,grid.AoA,grid.Lift,'spline','linear');
auxdata.interp.DragGridded = griddedInterpolant(grid.M,grid.AoA,grid.Drag,'spline','linear');
auxdata.interp.MomentGridded = griddedInterpolant(grid.M,grid.AoA,grid.Moment,'spline','linear');

%% Create Viscous Drag Interpolant

Aero_Visc1 = dlmread('VC3D_C_v_SPARTANFALCON_descent.dat');


interp.Cl_scattered_Viscousaero_1 = scatteredInterpolant(Aero_Visc1(:,1),Aero_Visc1(:,2),Aero_Visc1(:,3),Aero_Visc1(:,4));
interp.Cd_scattered_Viscousaero_1 = scatteredInterpolant(Aero_Visc1(:,1),Aero_Visc1(:,2),Aero_Visc1(:,3),Aero_Visc1(:,5));

MList_Visc1 = [0; unique(Aero_Visc1(:,1))]; % 0 included to prevent bad extrapolation

AoAList_Visc1 = unique(Aero_Visc1(:,2));

altList_Visc1 = unique(Aero_Visc1(:,3)); 

[Mgrid_Visc1,AOAgrid_Visc1,altgrid_Visc1] = ndgrid(MList_Visc1,AoAList_Visc1,altList_Visc1);

for i = 1:numel(Mgrid_Visc1)
    M_temp = Mgrid_Visc1(i);
    AoA_temp = AOAgrid_Visc1(i);
    alt_temp = altgrid_Visc1(i);
    
    I = cell(1, ndims(Mgrid_Visc1)); 
    [I{:}] = ind2sub(size(Mgrid_Visc1),i);
    
    Cl_temp_Viscous1 = interp.Cl_scattered_Viscousaero_1(M_temp,AoA_temp,alt_temp);
    Cd_temp_Viscous1= interp.Cd_scattered_Viscousaero_1(M_temp,AoA_temp,alt_temp);
    
    Cl_Grid_Visc1(I{(1)},I{(2)},I{(3)}) = Cl_temp_Viscous1;
    Cd_Grid_Visc1(I{(1)},I{(2)},I{(3)}) = Cd_temp_Viscous1;
    
end

auxdata.interp.Cl_gridded_Visc1 = griddedInterpolant(Mgrid_Visc1,AOAgrid_Visc1,altgrid_Visc1,Cl_Grid_Visc1,'linear','linear'); % Note: spline interp here was found to produce some odd results in in-between values
auxdata.interp.Cd_gridded_Visc1 = griddedInterpolant(Mgrid_Visc1,AOAgrid_Visc1,altgrid_Visc1,Cd_Grid_Visc1,'linear','linear');


%% Calculate Masses

mTotal = mSpartan + mRocket;
mEmpty = mRocket-mFuel;  %(kg)  %mass of the rocket (without fuel)


%% Assign Pitchover Conditions

%Define initial conditions at pitchover, these are assumed
h0 = 30;  
v0 = 15;    

gamma0 = deg2rad(89.9);    % set pitchover amount (start flight angle). This pitchover is 'free' movement, and should be kept small. 

mF = mEmpty+mSpartan;  %Assume that we use all of the fuel


alpha0 = 0; %Set initial angle of attack to 0

hMin1 = 1;   %Cannot go through the earth
hMax1 = 30000;  

vMin1 = 0; 
vMax1 = 3000;  

mMin1 = mEmpty;
mMax1 = mTotal;

if mode == 0
phiMin1 = -1.5;
phiMax1 = -0.5;
else
phiMin1 = -0.5;
phiMax1 = 0.5;
end
% phiMin1 = -0.5;
% phiMax1 = 0.5;

% zetaMin1 = -2*pi;
% zetaMax1 = 2*pi;

if mode == 0
zetaMin1 = -2*pi;
zetaMax1 = -pi/2; % constrained to launch over sea, to the west   
else
zetaMin1 = -2*pi;
zetaMax1 = 2*pi;
end

alphaMin1 = -deg2rad(5);
% alphaMax1 = deg2rad(2);
alphaMax1 = deg2rad(0);

ThrottleMin1 = 0.7;
ThrottleMax1 = 1;

dalphadt2Min1 = -0.1;
dalphadt2Max1 = 0.1;

lonMin = 2;         lonMax = 3;

gammaMin1 = deg2rad(-.1);
gammaMax1 = gamma0;
% This sets the control limits, this is second derivative of AoA
uMin1 = [-.0005, -0.05]; % AoAdot2 and thrustdot
uMax1 = [.0005, 0.05];

%-------------------------------------------
% Set up the problem bounds in SCALED units
%-------------------------------------------

tfMax1 	    = 300;     % large upper bound; do not choose Inf
		 
%-------------------------------------------------------------------------%
%---------- Provide Bounds and Guess in Each Phase of Problem ------------%
%-------------------------------------------------------------------------%

bounds.phase(1).initialtime.lower = 0;
bounds.phase(1).initialtime.upper = 0;

bounds.phase(1).finaltime.lower = 0;
bounds.phase(1).finaltime.upper = tfMax1;

bounds.phase(1).initialstate.lower = [h0, v0,  mF-1, gamma0, alpha0,  zetaMin1, dalphadt2Min1, lat0, lon0 ,ThrottleMin1];
bounds.phase(1).initialstate.upper = [h0, v0, mMax1, gamma0, alpha0, zetaMax1, dalphadt2Max1, lat0, lon0, ThrottleMax1 ];

bounds.phase(1).state.lower = [hMin1, vMin1, mF-1, gammaMin1, alphaMin1, zetaMin1, dalphadt2Min1, phiMin1, lonMin ,ThrottleMin1];
bounds.phase(1).state.upper = [ hMax1,  vMax1, mMax1, gammaMax1, alphaMax1, zetaMax1, dalphadt2Max1, phiMax1, lonMax, ThrottleMax1 ];

bounds.phase(1).finalstate.lower = [hMin1, vMin1, mF-1, gammaMin1, alphaMin1, zetaMin1, dalphadt2Min1, phiMin1, lonMin ,ThrottleMin1];
bounds.phase(1).finalstate.upper = [ hMax1,  vMax1, mMax1, gammaMax1, alphaMax1, zetaMax1, dalphadt2Max1, phiMax1, lonMax, ThrottleMax1 ];

bounds.phase(1).control.lower = uMin1;
bounds.phase(1).control.upper = uMax1;

bounds.phase(1).path.lower = 0;
bounds.phase(1).path.upper = 50000;

% Tie stages together
bounds.eventgroup(1).lower = [zeros(1,8)];
bounds.eventgroup(1).upper = [zeros(1,8)]; 

guess.phase(1).time =  [0; tfMax1];

guess.phase(1).state(:,1) = [h0; h0];
guess.phase(1).state(:,2) = [v0; 1500];
guess.phase(1).state(:,3) = [mMax1; mF];
guess.phase(1).state(:,4) = [gamma0; 0];
guess.phase(1).state(:,5) = [alpha0; 0];

if mode == 0
guess.phase(1).state(:,6) = [-pi;-pi];
else
guess.phase(1).state(:,6) = [0; 0];
end
guess.phase(1).state(:,7) = [0; 0];
guess.phase(1).state(:,8) = [lat0; lat0];
guess.phase(1).state(:,9) = [lon0; lon0];
guess.phase(1).state(:,10) = [1; 1];

guess.phase(1).control(:,1) = [0; 0];
guess.phase(1).control(:,2) = [0; 0];




%% Atmosphere Data %%======================================================
% Fetch atmospheric data and compute interpolation splines.

Atmosphere = dlmread('atmosphere.txt');
interp.Atmosphere = Atmosphere;
auxdata.interp.Atmosphere = interp.Atmosphere;

auxdata.interp.c_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,5)); % Calculate speed of sound using atmospheric data
auxdata.interp.rho_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,4)); % Calculate density using atmospheric data
auxdata.interp.p_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,3)); % Calculate density using atmospheric data
auxdata.interp.T0_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,2)); 
auxdata.interp.P0_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,3)); 

%% Import Vehicle and trajectory Config Data %%============================

run VehicleConfig.m
run TrajectoryConfig50kPa.m

auxdata.Stage3 = Stage3;
auxdata.Stage2 = Stage2;

%%
auxdata.Re   = 6371203.92;                     % Equatorial Radius of Earth (m)
auxdata.A = 62.77; %m^2

%% Third Stage Aerodynamic Data

auxdata.Aero3 = dlmread('ThirdStageAeroCoeffs.txt');

Aero3 = auxdata.Aero3;
auxdata.interp.Drag_interp3 = scatteredInterpolant(Aero3(:,1),Aero3(:,2),Aero3(:,5));
% 
auxdata.interp.Lift_interp3 = scatteredInterpolant(Aero3(:,1),Aero3(:,2),Aero3(:,6));

auxdata.interp.CP_interp3 = scatteredInterpolant(Aero3(:,1),Aero3(:,2),Aero3(:,7));

auxdata.interp.CN_interp3 = scatteredInterpolant(Aero3(:,1),Aero3(:,2),Aero3(:,4));

auxdata.interp.Max_AoA_interp3 = scatteredInterpolant(Aero3(:,1),Aero3(:,4),Aero3(:,2));


%% Conical Shock Data %%===================================================
% Import conical shock data and create interpolation splines 
shockdata = dlmread('ShockMat');
[MList_EngineOn,AOAList_EngineOn] = ndgrid(unique(shockdata(:,1)),unique(shockdata(:,2)));
M1_Grid = reshape(shockdata(:,3),[length(unique(shockdata(:,2))),length(unique(shockdata(:,1)))]).';
pres_Grid = reshape(shockdata(:,4),[length(unique(shockdata(:,2))),length(unique(shockdata(:,1)))]).';
temp_Grid = reshape(shockdata(:,5),[length(unique(shockdata(:,2))),length(unique(shockdata(:,1)))]).';
auxdata.interp.M1gridded = griddedInterpolant(MList_EngineOn,AOAList_EngineOn,M1_Grid,'spline','linear');
auxdata.interp.presgridded = griddedInterpolant(MList_EngineOn,AOAList_EngineOn,pres_Grid,'spline','linear');
auxdata.interp.tempgridded = griddedInterpolant(MList_EngineOn,AOAList_EngineOn,temp_Grid,'spline','linear');

%% Equivalence Ratio %%==========================================================
% Import engine data
auxdata.interp.engine_data = dlmread('ENGINEDATASORTED.txt');  % reads four columns; Mach no after conical shock, temp after conical shock, Isp, max equivalence ratio
engine_data = auxdata.interp.engine_data;

% Create uniform grid of Mach no. and temperature values. 
M_englist = unique(sort(engine_data(:,1))); % create unique list of Mach numbers from engine data
M_eng_interp = unique(sort(engine_data(:,1)));

T_englist = unique(sort(engine_data(:,2))); % create unique list of angle of attack numbers from engine data
T_eng_interp = unique(sort(engine_data(:,2)));

[grid.Mgrid_eng,grid.T_eng] =  ndgrid(M_eng_interp,T_eng_interp);

% Set the equivalence ratio interpolation region %-------------------------
% VERY IMPORTANT

% The interpolators have trouble with equivalence ratio because its equal
% to 1 over a certain Mach no. (causes error in interpolator, as the
% interpolator will find values of equivalence ratio < 1 where they should
% not exist)

% This makes anything outside of the region where it is actually changing
% extrapolate to over 1 (which is then set to 1 by RESTM12int)

% the the maximum of this to around where equivalence ratio stops changing,
% and check the end results

eq_data = [];
j=1;
for i = 1: length(engine_data(:,1))
    if engine_data(i,1) < 5.
        eq_data(j,:) = engine_data(i,:);
        j=j+1;
    end
end

auxdata.interp.equivalence = scatteredInterpolant(eq_data(:,1),eq_data(:,2),eq_data(:,4), 'linear');
grid.eq_eng = auxdata.interp.equivalence(grid.Mgrid_eng,grid.T_eng);
auxdata.interp.eqGridded = griddedInterpolant(grid.Mgrid_eng,grid.T_eng,grid.eq_eng,'linear','linear');

%% Isp data %-----------------------------------------

% gridIsp_eng is the spline interpolated data set created by
% engineint.m and engineinterpolator.exe
% 
% load gridIsp_eng
% grid.Isp_eng = gridIsp_eng;
% 
% % gridIsp_eng may have sections at which the Isp is 0. The following finds
% % these, and fills them in with linearly intepolated values.
% Isp_interpolator = scatteredInterpolant(engine_data(:,1),engine_data(:,2),engine_data(:,3));
% auxdata.interp.Isp_interpolator = Isp_interpolator;
% for i = 1:30 % must match engineint.m
%     for j= 1:30
%         % grid.Isp_eng(i,j) = polyvaln(p,[grid.Mgrid_eng(i,j) grid.T_eng(i,j)]);
%         if any(grid.Isp_eng(i,j)) == false
%             grid.Isp_eng(i,j) = Isp_interpolator(grid.Mgrid_eng(i,j), grid.T_eng(i,j)); % Linearly extrapolate for any data point which the Bivar spline could not solve
%         end
%     end
% end
% 
% auxdata.interp.IspGridded = griddedInterpolant(grid.Mgrid_eng,grid.T_eng,grid.Isp_eng,'spline','spline');
% auxdata.interp.IspGridded = griddedInterpolant(grid.Mgrid_eng,grid.T_eng,grid.Isp_eng,'linear','linear');


% auxdata.interp.Isprbf=rbfcreate([engine_data(:,1)'; engine_data(:,2)'], engine_data(:,3)' ,'RBFFunction', 'gaussian','RBFSmooth', 1);

% auxdata.interp.Isppoly = polyfitn([engine_data(:,1),engine_data(:,2)],engine_data(:,3),4);


%% Interpolate for Isp
% Note: no longer using Isp_grid, or files to generate that

% The engine data set is arranged to as to be difficult to interpolate
% effectively. Small imperfections in the interpolation scheme can have
% large effects on the optimal trajectory. To account for this, each 'line'
% of the C-rest data set is normalised so that the whole data set forms a
% regular 6x5 square, with even coordinate spacing. This is interpolated
% between using cubic interpolation. When this is interpolated for using a
% physical query point, the position of the data point in normalised
% coordinates is found by first finding the position of the query point
% relative to the existing data set. This position is converted to normalised 
% coordinates using the relative distance to the closest four data points. 

% Create coordinate system normalised to each 'line' in the engine data set
coords_noextrap = [];
temp = 1;

% create coordinates list for existing data set
for i = 2:7 % need to preallocate with extrapolation in mind
    for j = 2:6
        coords_noextrap(temp,1) = i;
        coords_noextrap(temp,2) = j;
        temp = temp+1;
        
        
    end
end

scattered_Isp_norm = scatteredInterpolant(coords_noextrap(:,1),coords_noextrap(:,2),engine_data(:,3));
scattered_M_norm = scatteredInterpolant(coords_noextrap(:,1),coords_noextrap(:,2),engine_data(:,1));
scattered_T_norm = scatteredInterpolant(coords_noextrap(:,1),coords_noextrap(:,2),engine_data(:,2));



% extrapolate in normalised coordinates and create grid
coords = [];
temp = 1;

for i = 1:8 % extrapolate by 1 'line'
    for j = 1:9
        coords(temp,1) = i;
        coords(temp,2) = j;
%         grid.Ispnorm(i,j) = engine_data(temp,3); % form an Isp grid at each coordinate
        grid.Ispnorm(i,j) = scattered_Isp_norm(i,j); % form an Isp grid at each coordinate
        M_List_data(temp) = scattered_M_norm(i,j);
        T_List_data(temp) = scattered_T_norm(i,j);
        temp = temp+1;

    end
end

M_List_data = M_List_data';
T_List_data = T_List_data';


[grid.normx,grid.normy] =  ndgrid(unique(coords(:,1)),unique(coords(:,2))); % form grids of normalised coordinates


% Create cubic spline on normalised coordinates
auxdata.interp.Ispnorm = griddedInterpolant(grid.normx,grid.normy,grid.Ispnorm,'cubic','linear'); % Cubic used as it is the best representation of the data

% Create the physical space over which to interpolate
M_eng_interp2 = min(engine_data(:,1)):0.1:max(engine_data(:,1));
% T_eng_interp2 = min(engine_data(:,2)):1:max(engine_data(:,2));
T_eng_interp2 = 220:1:max(engine_data(:,2));

[grid.Mgrid_eng,grid.T_eng] =  ndgrid(M_eng_interp2,T_eng_interp2);

size_grid = size(grid.Mgrid_eng);

Isp_grid= [];
M_list=[];
T_list=[];
Isp_norm_list=[];
temp = 1;

% Interpolate for M and T grids, this is done external of the main vehicle
% simulation for efficiency
for i = 1:size_grid(1)
    for j = 1:size_grid(2)
M_temp = grid.Mgrid_eng(i,j);
T_temp = grid.T_eng(i,j);

in_dex = [0 0 0 0];

% Determine if the qeury point is within any four points of the data set,
% and if so, which four
for polyindex_x = 1:7
    for polyindex_y = 1:8
        % Create polygon fof four data points
        cells_index = [polyindex_y+9*(polyindex_x-1)  polyindex_y+1+9*(polyindex_x-1)  polyindex_y+1+9*(polyindex_x) polyindex_y+9*(polyindex_x)];
        % Search
        in = inpolygon(M_temp,T_temp,M_List_data(cells_index),T_List_data(cells_index));
        if in
            in_dex = cells_index; % record cells index for polygon if query point inside
        end
    end
end

if in_dex == [0 0 0 0]
Isp_grid(i,j) = 0; % if not inside data set set point to 0

else
    % if inside a polygon of data points
    % Define a polygon of the data points
pt0 = [M_temp T_temp 0]; %query pt
pt1 = [M_List_data(in_dex(1)) T_List_data(in_dex(1)) 0]; %data pt 1
pt2 = [M_List_data(in_dex(2)) T_List_data(in_dex(2))  0];
pt3 = [M_List_data(in_dex(3)) T_List_data(in_dex(3))  0];
pt4 = [M_List_data(in_dex(4)) T_List_data(in_dex(4)) 0];

% Calculate minimum distances from query pt to a line between each 'side'
% of data pts
      a1 = pt1 - pt2;
      b1 = pt0 - pt2;
      d1 = norm(cross(a1,b1)) / norm(a1);

      a2 = pt3 - pt4;
      b2 = pt0 - pt4;
      d2 = norm(cross(a2,b2)) / norm(a2);
      
      a3 = pt1 - pt4;
      b3 = pt0 - pt4;
      d3 = norm(cross(a3,b3)) / norm(a3);
      
      a4 = pt2 - pt3;
      b4 = pt0 - pt3;
      d4 = norm(cross(a4,b4)) / norm(a4);
      
      % Calculate normalised coordinate as the distance fraction from the southwest data pt to the query pt in the x
      % and y directions, added to the coordinate of the southwest data pt
x_norm(i,j) = coords(in_dex(1),1) + d1/(d1+d2);
y_norm(i,j) = coords(in_dex(1),2) + d3/(d3+d4);  

Isp_grid(i,j) = auxdata.interp.Ispnorm(x_norm(i,j),y_norm(i,j)); % interpolate for Isp and create grid corresponding to M and T

% Define lists for extrapolation purposes
M_list(temp) = M_temp;
T_list(temp) = T_temp;
Isp_norm_list(temp) = Isp_grid(i,j);
temp = temp+1;
end

    end
end

% Use the interpolated data set to generate a scattered interpolant for
% extrapolation
auxdata.interp.Ispscatterednorm = scatteredInterpolant(M_list',T_list',Isp_norm_list','linear','nearest');


% % Extrapolate on any point outside of data set
for i = 1:size_grid(1)
    for j = 1:size_grid(2)
        if Isp_grid(i,j) == 0
            
            Isp_grid(i,j) = auxdata.interp.Ispscatterednorm(grid.Mgrid_eng(i,j),grid.T_eng(i,j));
        end
    end
end
% contourf(grid.Mgrid_eng,grid.T_eng,Isp_grid,100)
% Create spline for use in vehicle model.
auxdata.interp.IspGridded = griddedInterpolant(grid.Mgrid_eng,grid.T_eng,Isp_grid,'spline','spline');
%% Aerodynamic Data

% Fetch aerodynamic data and compute interpolation splines.
% Calculate the flap deflection necessary for trim.
% Each set of aero corresponds to a different CG. 



Viscousaero_EngineOn = importdata('VC3D_viscousCoefficients_ascent.dat');
Viscousaero_EngineOff = importdata('VC3D_viscousCoefficients_descent.dat');

% These aerodynamic datasets have been created in ClicCalcCGVar.m

T_L = -1.327; % Thrust location in z, average (m), measured from CREO

% Full of fuel, with third stage

CG_z = 0.049; % Centre of gravity location in z, calculated from Creo

aero1.aero_EngineOff = importdata('SPARTANaero15.228');
aero1.flapaero = importdata('SPARTANaeroFlaps15.228');
aero1.aero_EngineOn = importdata('SPARTANaeroEngineOn15.228');
aero1.aero_Engine = importdata('SPARTANEngine15.228');
aero1.Viscousaero_EngineOff = Viscousaero_EngineOff;
aero1.Viscousaero_EngineOn = Viscousaero_EngineOn;

[auxdata.interp.Cl_spline_EngineOff.fullFuel,auxdata.interp.Cd_spline_EngineOff.fullFuel,...
    auxdata.interp.Cl_spline_EngineOn.fullFuel,auxdata.interp.Cd_spline_EngineOn.fullFuel,...
    auxdata.interp.flap_spline_EngineOff.fullFuel,auxdata.interp.flap_spline_EngineOn.fullFuel,...
    auxdata.T_spline_Rear,auxdata.Fd_spline_NoEngine,auxdata.Cd_spline_ViscousEngineOff,auxdata.Cd_spline_ViscousEngineOn,...
    auxdata.L_spline_Rear,auxdata.T_spline] = AeroInt(aero1,auxdata,T_L,CG_z);


% Cylindrical fuel tanks depleted, with third stage
CG_z = 0.077;

aero2.aero_EngineOff = importdata('SPARTANaero15.1557');
aero2.flapaero = importdata('SPARTANaeroFlaps15.1557');
aero2.aero_EngineOn = importdata('SPARTANaeroEngineOn15.1557');
aero2.aero_Engine = importdata('SPARTANEngine15.1557');
aero2.Viscousaero_EngineOff = Viscousaero_EngineOff;
aero2.Viscousaero_EngineOn = Viscousaero_EngineOn;

[auxdata.interp.Cl_spline_EngineOff.cylTankEnd,auxdata.interp.Cd_spline_EngineOff.cylTankEnd,auxdata.interp.Cl_spline_EngineOn.cylTankEnd,auxdata.interp.Cd_spline_EngineOn.cylTankEnd,auxdata.interp.flap_spline_EngineOff.cylTankEnd,auxdata.interp.flap_spline_EngineOn.cylTankEnd] = AeroInt(aero2,auxdata,T_L,CG_z);



% Empty, with third stage. 


CG_z = 0.086;

aero3.aero_EngineOff = importdata('SPARTANaero15.727');
aero3.flapaero = importdata('SPARTANaeroFlaps15.727');
aero3.aero_EngineOn = importdata('SPARTANaeroEngineOn15.727');
aero3.aero_Engine = importdata('SPARTANEngine15.727');
aero3.Viscousaero_EngineOff = Viscousaero_EngineOff;
aero3.Viscousaero_EngineOn = Viscousaero_EngineOn;

[auxdata.interp.Cl_spline_EngineOff.noFuel,auxdata.interp.Cd_spline_EngineOff.noFuel,auxdata.interp.Cl_spline_EngineOn.noFuel,auxdata.interp.Cd_spline_EngineOn.noFuel,auxdata.interp.flap_spline_EngineOff.noFuel,auxdata.interp.flap_spline_EngineOn.noFuel] = AeroInt(aero3,auxdata,T_L,CG_z);

% Flyback, empty without third stage, empty.

CG_z = -0.2134; % calculated fom CREO

aero4.aero_EngineOff = importdata('SPARTANaero15.1255');
aero4.flapaero = importdata('SPARTANaeroFlaps15.1255');
aero4.aero_EngineOn = importdata('SPARTANaeroEngineOn15.1255');
aero4.aero_Engine = importdata('SPARTANEngine15.1255');
aero4.Viscousaero_EngineOff = Viscousaero_EngineOff;
aero4.Viscousaero_EngineOn = Viscousaero_EngineOn;

[auxdata.interp.Cl_spline_EngineOff.noThirdStageEmpty,...
    auxdata.interp.Cd_spline_EngineOff.noThirdStageEmpty,auxdata.interp.Cl_spline_EngineOn.noThirdStageEmpty,...
    auxdata.interp.Cd_spline_EngineOn.noThirdStageEmpty,auxdata.interp.flap_spline_EngineOff.noThirdStageEmpty,...
    auxdata.interp.flap_spline_EngineOn.noThirdStageEmpty] = AeroInt(aero4,auxdata,T_L,CG_z);

% Flyback, empty without third stage, conical fuel tank full.
CG_z = -0.183; % calculated fom CREO

aero5.aero_EngineOff = importdata('SPARTANaero14.297');
aero5.flapaero = importdata('SPARTANaeroFlaps14.297');
aero5.aero_EngineOn = importdata('SPARTANaeroEngineOn14.297');
aero5.aero_Engine = importdata('SPARTANEngine14.297');
aero5.Viscousaero_EngineOff = Viscousaero_EngineOff;
aero5.Viscousaero_EngineOn = Viscousaero_EngineOn;

[auxdata.interp.Cl_spline_EngineOff.noThirdStagecylTankEnd,...
    auxdata.interp.Cd_spline_EngineOff.noThirdStagecylTankEnd,auxdata.interp.Cl_spline_EngineOn.noThirdStagecylTankEnd,...
    auxdata.interp.Cd_spline_EngineOn.noThirdStagecylTankEnd,auxdata.interp.flap_spline_EngineOff.noThirdStagecylTankEnd,...
    auxdata.interp.flap_spline_EngineOn.noThirdStagecylTankEnd] = AeroInt(aero5,auxdata,T_L,CG_z);


clear aero1
clear aero2
clear aero3
clear aero4
clear aero5
clear Viscousaero_EngineOn
clear Viscousaero_EngineOff


%% Import Bounds %%========================================================

if mode == 0
latMin2 = -1.5;  latMax2 = -0.5;
else
latMin2 = -0.5;  latMax2 = 0.5;
end

% latMin2 = -0.5;  latMax2 = 0.5;



aoaMin21 = 0;  aoaMax21 = 10*pi/180;

if mode == 0
    bankMin21 = -90*pi/180; bankMax21 =   90*pi/180; 
else
    if returnMode == 0
    bankMin21 = -1*pi/180; bankMax21 =   1*pi/180;    
    else
    bankMin21 = -1*pi/180; bankMax21 =   90*pi/180;
    end
end

    
% if returnMode == 0
% bankMin21 = -1*pi/180; bankMax21 =   1*pi/180;    
% else
% bankMin21 = -1*pi/180; bankMax21 =   90*pi/180;
% end

if mode == 0
    zetaMin21 = -4*pi;
    zetaMax21 = 2*pi;
else
    zetaMin21 = Stage2.Bounds.zeta(1);
    zetaMax21 = Stage2.Bounds.zeta(2);
end

% Primal Bounds
bounds.phase(2).state.lower = [Stage2.Bounds.Alt(1), lonMin, latMin2, Stage2.Bounds.v(1), Stage2.Bounds.gamma(1), zetaMin21, aoaMin21, bankMin21, Stage2.Bounds.mFuel(1)];
bounds.phase(2).state.upper = [Stage2.Bounds.Alt(2), lonMax, latMax2, Stage2.Bounds.v(2), deg2rad(15), zetaMax21, aoaMax21, bankMax21, Stage2.Bounds.mFuel(2)];

% Initial States
if returnMode == 1
   bounds.phase(2).initialstate.lower = [Stage2.Bounds.Alt(1),lonMin, latMin2, Stage2.Bounds.v(1), Stage2.Bounds.gamma(1), zetaMin21, aoaMin21, 0, Stage2.Initial.mFuel] ;
bounds.phase(2).initialstate.upper = [Stage2.Bounds.Alt(2),lonMax, latMax2, Stage2.Bounds.v(2), deg2rad(15), zetaMax21, aoaMax21, 0, Stage2.Initial.mFuel];
 
else
bounds.phase(2).initialstate.lower = [Stage2.Bounds.Alt(1),lonMin, latMin2, Stage2.Bounds.v(1), Stage2.Bounds.gamma(1), zetaMin21, aoaMin21, bankMin21, Stage2.Initial.mFuel] ;
bounds.phase(2).initialstate.upper = [Stage2.Bounds.Alt(2),lonMax, latMax2, Stage2.Bounds.v(2), deg2rad(15), zetaMax21, aoaMax21, bankMax21, Stage2.Initial.mFuel];
end

% bounds.phase(2).initialstate.lower = [Stage2.Bounds.Alt(1),lon0, lat0, 1500, Stage2.Bounds.gamma(1), Stage2.Bounds.zeta(1), aoaMin21, bankMin21, Stage2.Initial.mFuel] ;
% bounds.phase(2).initialstate.upper = [Stage2.Bounds.Alt(2),lon0, lat0, 1500, deg2rad(15), Stage2.Bounds.zeta(2), aoaMax21, bankMax21, Stage2.Initial.mFuel];

% End States
% End bounds are set slightly differently, to encourage an optimal solution
bounds.phase(2).finalstate.lower = [20000, lonMin, latMin2, 2300, 0, zetaMin21, aoaMin21, 0, Stage2.End.mFuel];
if mode == 90
bounds.phase(2).finalstate.upper = [45000, lonMax, latMax2, Stage2.Bounds.v(2), deg2rad(0.5), zetaMax21, aoaMax21, 0, Stage2.Initial.mFuel];
else
bounds.phase(2).finalstate.upper = [45000, lonMax, latMax2, Stage2.Bounds.v(2), deg2rad(20), zetaMax21, aoaMax21, 0, Stage2.Initial.mFuel];
end
% bounds.phase(2).finalstate.lower = [34000, lonMin, latMin2, 2300, 0, Stage2.Bounds.zeta(1), aoaMin21, 0, Stage2.End.mFuel];
% bounds.phase(2).finalstate.upper = [45000, lonMax, latMax2, Stage2.Bounds.v(2), 0, Stage2.Bounds.zeta(2), aoaMax21, 0, Stage2.Initial.mFuel];
%  disp('gamma set to 0')
 
% Control Bounds
bounds.phase(2).control.lower = [deg2rad(-.5), deg2rad(-1)];
bounds.phase(2).control.upper = [deg2rad(.5), deg2rad(1)];

% Time Bounds
bounds.phase(2).initialtime.lower = 0;
bounds.phase(2).initialtime.upper = Stage2.Bounds.time(2);
bounds.phase(2).finaltime.lower = Stage2.Bounds.time(1);
bounds.phase(2).finaltime.upper = Stage2.Bounds.time(2);


%% Define Path Constraints
% Path bounds, defined in Continuous function.
% These limit the dynamic pressure.
% if mode == 1 || mode == 14 || mode == 15 || mode == 2

% if mode ==90
%         bounds.phase(2).path.lower = [49500];
%     bounds.phase(2).path.upper = [50500];
% else
    bounds.phase(2).path.lower = [0,0];
    bounds.phase(2).path.upper = [50000,8000];
% end
% elseif mode ==3 || mode == 32
%         bounds.phase(2).path.lower = [49970];
%     bounds.phase(2).path.upper = [50030];
% end
if mode == 90
bounds.phase(2).integral.lower = 0;
bounds.phase(2).integral.upper = 10000000;
end
%%  Guess =================================================================
% Set the initial guess. This can have a significant effect on the final
% solution, even for a well defined problem. 

% guess.phase(2).state(:,1)   = [24000;35000];
if mode == 0
guess.phase(2).state(:,2)   = [lon0;lon0];
guess.phase(2).state(:,3)   = [lat0;lat0-0.1]; 
else
guess.phase(2).state(:,2)   = [2.53;2.5368];
guess.phase(2).state(:,3)   = [-0.269;-0.10];
end

guess.phase(2).state(:,4)   = Stage2.Guess.v.';
guess.phase(2).state(:,5)   = Stage2.Guess.gamma.';

if mode == 0
    guess.phase(2).state(:,6)   = [-pi;-pi];
else
guess.phase(2).state(:,6)   = Stage2.Guess.zeta.';
end

guess.phase(2).state(:,7)   = [2*pi/180; 5*pi/180];
guess.phase(2).state(:,8)   = [deg2rad(10);deg2rad(10)];
guess.phase(2).state(:,9) 	= [Stage2.Initial.mFuel; 100];

guess.phase(2).control      = [[0;0],[0;0]];
guess.phase(2).time          = [0;650];

if mode == 90
guess.phase(2).integral = 0
end
% Tie stages together
bounds.eventgroup(2).lower = [zeros(1,9)];
bounds.eventgroup(2).upper = [zeros(1,9)]; 

%% Flyback
tfMin = 0;            tfMax = 5000;
altMin = 10;  altMax = 70000;
speedMin = 10;        speedMax = 5000;
fpaMin = -80*pi/180;  fpaMax =  80*pi/180;

if mode == 0
aziMin = -4*pi; aziMax =  4*pi;
else
aziMin = 60*pi/180; aziMax =  500*pi/180;
end

mFuelMin = 0; mFuelMax = 500;

if mode == 0
bankMin21 = -90*pi/180; bankMax21 =   90*pi/180;  
else
bankMin21 = -10*pi/180; bankMax21 =   90*pi/180;
end

throttleMin = 0; throttleMax = 1;

bounds.phase(3).initialtime.lower = 0;
bounds.phase(3).initialtime.upper = 3000;
bounds.phase(3).finaltime.lower = 400;
bounds.phase(3).finaltime.upper = 4000;

% Initial Bounds

 bounds.phase(3).initialstate.lower = [altMin, lonMin, latMin2, speedMin, fpaMin, aziMin, aoaMin21, 0, mFuelMin, throttleMin];
bounds.phase(3).initialstate.upper = [altMax, lonMax, latMax2, speedMax, fpaMax, aziMax, aoaMax21, 0, mFuelMax, throttleMax];   

% State Bounds
if returnMode == 0
bounds.phase(3).state.lower = [altMin, lonMin, latMin2, speedMin, fpaMin, aziMin, aoaMin21, bankMin21, mFuelMin, throttleMin];
bounds.phase(3).state.upper = [altMax, lonMax, latMax2, speedMax, fpaMax, aziMax, aoaMax21, bankMax21, 1, throttleMax];
else
bounds.phase(3).state.lower = [altMin, lonMin, latMin2, speedMin, fpaMin, aziMin, aoaMin21, bankMin21, mFuelMin, throttleMin];
bounds.phase(3).state.upper = [altMax, lonMax, latMax2, speedMax, fpaMax, aziMax, aoaMax21, bankMax21, mFuelMax, throttleMax];
end
% End State Bounds
if returnMode == 0
bounds.phase(3).finalstate.lower = [altMin, lonMin, latMin2, speedMin, fpaMin, aziMin, aoaMin21, bankMin21, mFuelMin, throttleMin];
bounds.phase(3).finalstate.upper = [altMax, lonMax, latMax2, speedMax, fpaMax, aziMax, aoaMax21, bankMax21, mFuelMax, throttleMax];
else
% bounds.phase(3).finalstate.lower = [altMin, lon0, lat0, speedMin, deg2rad(-20), aziMin, aoaMin21, bankMin21, Stage2.End.mFuel, throttleMin];
% bounds.phase(3).finalstate.upper = [altMax, lon0, lat0, speedMax, deg2rad(30), aziMax, aoaMax21, bankMax21, Stage2.End.mFuel, throttleMax];
bounds.phase(3).finalstate.lower = [altMin, lon0, lat0, speedMin, fpaMin, aziMin, aoaMin21, bankMin21, mFuelMin, throttleMin];
bounds.phase(3).finalstate.upper = [1000, lon0, lat0, speedMax, fpaMax, aziMax, aoaMax21, bankMax21, mFuelMax, throttleMax];
end
% Control Bounds
if returnMode == 1 % these being different has no purpose but to keep my ascent results consistent
 bounds.phase(3).control.lower = [deg2rad(-.5), deg2rad(-1), -.1];
bounds.phase(3).control.upper = [deg2rad(.5), deg2rad(1), .1];   
else
bounds.phase(3).control.lower = [deg2rad(-.5), deg2rad(-1), -.2];
bounds.phase(3).control.upper = [deg2rad(.5), deg2rad(1), .2];
end
% Path Bounds
bounds.phase(3).path.lower = [0,0];
bounds.phase(3).path.upper = [50000,8000];

bounds.eventgroup(3).lower = [0]; 
bounds.eventgroup(3).upper = [0]; 

% Guess
tGuess              = [440; 1500];
altGuess            = [35000; 100];

if mode == 0
lonGuess            = [lon0; lon0-.1*pi/180];
latGuess            = [lat0-0.1;lat0];
else
lonGuess            = [lon0; lon0-.1*pi/180];
latGuess            = [-0.11;-0.10-0*pi/180];
end

speedGuess          = [3000; 10];
fpaGuess            = [0; 0];

if mode == 0
aziGuess            = [-pi; 0]; 
else
aziGuess            = [deg2rad(97); deg2rad(270)];
end

aoaGuess            = [8*pi/180; 5*pi/180];
bankGuess           = [89*pi/180; 89*pi/180];
mFuelGuess          = [100; mFuelMin];
if returnMode == 1
guess.phase(3).state   = [altGuess, lonGuess, latGuess, speedGuess, fpaGuess, aziGuess, aoaGuess, bankGuess, mFuelGuess,[1.;1.]];
else
guess.phase(3).state   = [altGuess, lonGuess, latGuess, speedGuess, fpaGuess, aziGuess, aoaGuess, bankGuess, mFuelGuess,[0.;0.]];
end
guess.phase(3).control = [[0;0],[0;0],[0;0]];
guess.phase(3).time    = tGuess;


%% Third Stage

altMin3 = 30000;  altMax3 = 84000;
if mode == 0
phiMin3 = -1.5;         phiMax3 = -0.5;   
else
phiMin3 = -0.5;         phiMax3 = 0.5;
end
vMin3 = 10;        vMax3 = 8000;
gammaMin3=deg2rad(-5);  gammaMax3 =  deg2rad(30);

if mode == 0
zetaMin3 = deg2rad(-180); zetaMax3 =  deg2rad(0); 
else
zetaMin3 = deg2rad(80); zetaMax3 =  deg2rad(120);
end

aoaMin3 = 0;  aoaMax3= deg2rad(20);
aoadotMin3 = -deg2rad(1);
aoadotMax3 = deg2rad(1);

%-------------------------------------------------------------------%
%--------------- Set Up Problem Using Data Provided Above ----------%
%-------------------------------------------------------------------%
bounds.phase(4).initialtime.lower = 0;
bounds.phase(4).initialtime.upper = 10000;
bounds.phase(4).finaltime.lower = 1;
bounds.phase(4).finaltime.upper = 10000;


bounds.phase(4).initialstate.lower = [altMin3,vMin3, 0,  auxdata.Stage3.mTot, aoaMin3, phiMin3, zetaMin3];
bounds.phase(4).initialstate.upper = [altMax3, vMax3, gammaMax3, auxdata.Stage3.mTot, aoaMax3, phiMax3, zetaMax3];

bounds.phase(4).state.lower = [altMin3,vMin3, gammaMin3, 0, aoaMin3, phiMin3, zetaMin3];
bounds.phase(4).state.upper = [altMax3, vMax3, gammaMax3, auxdata.Stage3.mTot, aoaMax3, phiMax3, zetaMax3];

bounds.phase(4).finalstate.lower = [altMin3, vMin3, 0, 0, 0, phiMin3, zetaMin3];
bounds.phase(4).finalstate.upper = [altMax3, vMax3, gammaMax3, auxdata.Stage3.mTot, 0, phiMax3, zetaMax3];

bounds.phase(4).control.lower = [aoadotMin3];
bounds.phase(4).control.upper = [aoadotMax3];

bounds.phase(4).path.lower = [-deg2rad(8), -inf];
bounds.phase(4).path.upper = [deg2rad(8), 0];


bounds.eventgroup(4).lower = [zeros(1,7) 90000 0];
bounds.eventgroup(4).upper = [zeros(1,6) 1000 566000 0];

tGuess              = [0; 150];
altGuess            = [35000; 60000];
vGuess          = [2700; 6000];
gammaGuess            = [0; deg2rad(10)];
mGuess              = [3300; 2000];
aoaGuess            = [deg2rad(20); deg2rad(20)];



if mode == 0
phiGuess = [lat0-0.1;lat0-0.1];
zetaGuess = [-deg2rad(97);-deg2rad(97)];  
else
phiGuess = [-0.11;-0.11];
zetaGuess = [deg2rad(97);deg2rad(97)];
end

guess.phase(4).state   = [altGuess, vGuess, gammaGuess, mGuess, aoaGuess, phiGuess, zetaGuess];
guess.phase(4).control = [0;0];
guess.phase(4).time    = tGuess;

%%
%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method       = 'hp-LiuRao-Legendre'; 
% mesh.method       = 'hp-LiuRao'; 
%  mesh.method       = 'hp-DarbyRao';
% mesh.maxiterations = 10;
if mode == 90
  mesh.maxiterations = 6;  
elseif returnMode == 1
 mesh.maxiterations = 4;   
else
mesh.maxiterations = 3;
end
mesh.colpointsmin = 8;
mesh.colpointsmax = 500;
mesh.tolerance    = 1e-5;

%--------------------------------------------------------%
%---------- Configure Setup Using the information provided ---------%
%-------------------------------------------------------------------%

setup.name                           = 'SPARTAN-Combined';
setup.functions.continuous           = @CombinedContinuous;
setup.functions.endpoint             = @CombinedEndpoint;
setup.auxdata                        = auxdata;
setup.bounds                         = bounds;
setup.guess                          = guess;
setup.mesh                           = mesh;
setup.displaylevel                   = 2;
setup.nlp.solver                     = 'ipopt';
setup.nlp.ipoptoptions.linear_solver = 'ma57';

setup.nlp.ipoptoptions.maxiterations = 150;

setup.derivatives.supplier           = 'sparseFD';

setup.derivatives.derivativelevel    = 'first';
% setup.scales.method                  = 'automatic-bounds';
setup.method                         = 'RPM-Differentiation';
setup.scales.method                  = 'automatic-guessUpdate';
setup.derivatives.dependencies      = 'full';



%-------------------------------------------------------------------%
%------------------- Solve Problem Using GPOPS2 --------------------%
%-------------------------------------------------------------------%

%% Parallel Loop

num_it = 4; % Define number of parallel iterations

% Create variable setup structure
%% for interaction testing
if mode == 99
    modes_list = [2 3 4]; % define modes to generate interactions for
    runs1 = [-1 -1  zeros(1,length(modes_list))]; % defines the number of runs to perform. Need to define two -1 and two 1
    runs2 = [-1  1 zeros(1,length(modes_list))]; % defines the number of runs to perform. Need to define two -1 and two 1
    runs3 = [ 1 1 zeros(1,length(modes_list))]; % defines the number of runs to perform. Need to define two -1 and two 1
    
    % calculate all possible permutations where 2 of the given modes are
    % chosen
    k = length(modes_list);
    nk1 = nchoosek(runs1,k);
    p1=zeros(0,k);
    for i=1:size(nk1,1),
        pi1 = perms(nk1(i,:));
        p1 = unique([p1; pi1],'rows');
    end
    % calculate all possible permutations where 2 of the given modes are
    % chosen
    k = length(modes_list);
    nk2 = nchoosek(runs2,k);
    p2=zeros(0,k);
    for i=1:size(nk2,1),
        pi2 = perms(nk2(i,:));
        p2 = unique([p2; pi2],'rows');
    end
    % calculate all possible permutations where 2 of the given modes are
    % chosen
    k = length(modes_list);
    nk3 = nchoosek(runs3,k);
    p3=zeros(0,k);
    for i=1:size(nk3,1),
        pi3 = perms(nk3(i,:));
        p3 = unique([p3; pi3],'rows');
    end
    p = [p1;p2;p3];
   
     auxdata.p = p;
% Step through modes
    for j = 1:size(p,1)
        namelist{j} = num2str(j);
        setup_variations{j} = setup;
        if p(j,1) == 1
                
                setup_variations{j}.bounds.phase(1).path.upper = q_vars(end-1);
                setup_variations{j}.bounds.phase(2).path.upper = q_vars(end-1);
                setup_variations{j}.bounds.phase(3).path.upper = q_vars(end-1);
        elseif p(j,1) == -1
 
                setup_variations{j}.bounds.phase(1).path.upper = q_vars(2);
                setup_variations{j}.bounds.phase(2).path.upper = q_vars(2);
                setup_variations{j}.bounds.phase(3).path.upper = q_vars(2);
        end
        if p(j,2) == 1
          
                setup_variations{j}.auxdata.Ispmod = Isp_vars(end-1);
        elseif p(j,2) == -1
              
                setup_variations{j}.auxdata.Ispmod = Isp_vars(2);
        end
        if p(j,3) == 1
           
                setup_variations{j}.auxdata.Cdmod = Cd_vars(end-1);
        elseif p(j,3) == -1
             
                setup_variations{j}.auxdata.Cdmod = Cd_vars(2);
        end
%         if mode == 5 % not yet implemented
%                 setup_variations{j} = setup;
%                 setup_variations{j}.auxdata.vCdmod = vCd_vars(end);
%         end
%         if mode == 8 %
%                 setup_variations{j} = setup;
%                 setup_variations{j}.bounds.phase(4).initialstate.lower(4) = auxdata.Stage3.mTot*m3_vars(end);
%                 setup_variations{j}.bounds.phase(4).initialstate.upper(4) = auxdata.Stage3.mTot*m3_vars(end);
%                 setup_variations{j}.bounds.phase(4).state.upper(4) = auxdata.Stage3.mTot*m3_vars(end);
%                 setup_variations{j}.bounds.phase(4).finalstate.upper(4) = auxdata.Stage3.mTot*m3_vars(end);
%         end
%         if mode == 9 %
%                 setup_variations{j} = setup;
%                 setup_variations{j}.auxdata.Isp3mod = Isp3_vars(end);
%         end
%         if mode == 10 %
%                 setup_variations{j} = setup;
%                 setup_variations{j}.auxdata.Stage2.mStruct = auxdata.Stage2.mStruct*mSPARTAN_vars(end);
%         end
%         if mode == 11 %
%                 setup_variations{j} = setup;
%                 setup_variations{j}.bounds.phase(2).initialstate.lower(9) = Stage2.Initial.mFuel*mFuel_vars(end);
%                 setup_variations{j}.bounds.phase(2).initialstate.upper(9) = Stage2.Initial.mFuel*mFuel_vars(end);
%                 setup_variations{j}.bounds.phase(2).state.upper(9) = Stage2.Initial.mFuel*mFuel_vars(end);
%                 setup_variations{j}.bounds.phase(2).finalstate.upper(9) = Stage2.Initial.mFuel*mFuel_vars(end);
% 
%         end

    end
end

%% for independent variable testing
if mode == 1 || mode == 0
    setup_variations{1} = setup;
elseif mode == 90
    setup_variations{1} = setup;
elseif mode == 2
    for i = 1:length(q_vars)  
        setup_variations{i} = setup;
        setup_variations{i}.bounds.phase(1).path.upper = q_vars(i);
        setup_variations{i}.bounds.phase(2).path.upper = q_vars(i);
        setup_variations{i}.bounds.phase(3).path.upper = q_vars(i);
    end
%     setup_variations{2}.guess.phase(2).state(:,4) = [1500;2800]; % modify guess for badly converging solution
elseif mode == 3
    for i = 1:length(Isp_vars)  
        setup_variations{i} = setup;
        setup_variations{i}.auxdata.Ispmod = Isp_vars(i);
    end
%     setup_variations{2}.guess.phase(2).state(:,4) = [1500;2800]; % modify guess for badly converging solution
%     if returnMode == 1
%         setup_variations{1}.guess.phase(2).state(:,4) = [1500;2700]; % modify guess for badly converging solution
%         setup_variations{2}.guess.phase(2).state(:,4) = [1500;2700]; % modify guess for badly converging solution
%         setup_variations{4}.guess.phase(2).state(:,4) = [1500;2700]; % modify guess for badly converging solution
%         setup_variations{5}.guess.phase(2).state(:,4) = [1500;2700]; % modify guess for badly converging solution
%     end
elseif mode == 4
    for i = 1:length(Cd_vars)  
        setup_variations{i} = setup;
        setup_variations{i}.auxdata.Cdmod = Cd_vars(i);
    end
    
elseif mode == 5 % not yet implemented
    for i = 1:length(vCd_vars)  
        setup_variations{i} = setup;
        setup_variations{i}.auxdata.vCdmod = vCd_vars(i);
    end
    elseif mode == 6 % not yet implemented
    for i = 1:length(FirstStagem_vars)  
        setup_variations{i} = setup;
        
%         setup_variations{i}.bounds.phase(1).initialstate.upper(3) = mMax1+mRocket*(FirstStagem_vars(i)-1);
%         setup_variations{i}.bounds.phase(1).state.upper(3) = mMax1+mRocket*(FirstStagem_vars(i)-1);
        setup_variations{i}.bounds.phase(1).finalstate.lower(3) = mF + mEmpty*(FirstStagem_vars(i)-1);
%         setup_variations{i}.bounds.phase(1).finalstate.upper(3) = mMax1+mRocket*(FirstStagem_vars(i)-1);
    end
    elseif mode == 7 %
    for i = 1:length(Cd3_vars)  
        setup_variations{i} = setup;
        setup_variations{i}.auxdata.Cd3mod = Cd3_vars(i);
    end
elseif mode == 8 %
    for i = 1:length(m3_vars)  
        setup_variations{i} = setup;
        setup_variations{i}.bounds.phase(4).initialstate.lower(4) = auxdata.Stage3.mTot*m3_vars(i);
        setup_variations{i}.bounds.phase(4).initialstate.upper(4) = auxdata.Stage3.mTot*m3_vars(i);
        setup_variations{i}.bounds.phase(4).state.upper(4) = auxdata.Stage3.mTot*m3_vars(i);
        setup_variations{i}.bounds.phase(4).finalstate.upper(4) = auxdata.Stage3.mTot*m3_vars(i);
        
        setup_variations{i}.bounds.phase(1).initialstate.lower(3) = setup.bounds.phase(1).initialstate.lower(3) + auxdata.Stage3.mTot*(m3_vars(i)-1);
        setup_variations{i}.bounds.phase(1).initialstate.upper(3) = setup.bounds.phase(1).initialstate.upper(3) + auxdata.Stage3.mTot*(m3_vars(i)-1);
        setup_variations{i}.bounds.phase(1).state.lower(3) = setup.bounds.phase(1).state.lower(3) + auxdata.Stage3.mTot*(m3_vars(i)-1);
        setup_variations{i}.bounds.phase(1).state.upper(3) = setup.bounds.phase(1).state.upper(3) + auxdata.Stage3.mTot*(m3_vars(i)-1);
        setup_variations{i}.bounds.phase(1).finalstate.lower(3) = setup.bounds.phase(1).finalstate.lower(3) + auxdata.Stage3.mTot*(m3_vars(i)-1);
        setup_variations{i}.bounds.phase(1).finalstate.upper(3) = setup.bounds.phase(1).finalstate.upper(3) + auxdata.Stage3.mTot*(m3_vars(i)-1);
    
        setup_variations{i}.auxdata.Stage3.mTot = auxdata.Stage3.mTot*m3_vars(i);
    end
elseif mode == 9 %
    for i = 1:length(Isp3_vars)  
        setup_variations{i} = setup;
        setup_variations{i}.auxdata.Isp3mod = Isp3_vars(i);
    end
elseif mode == 10 %
    for i = 1:length(mSPARTAN_vars)  
        setup_variations{i} = setup;
        setup_variations{i}.auxdata.Stage2.mStruct = auxdata.Stage2.mStruct*mSPARTAN_vars(i);
        setup_variations{i}.bounds.phase(1).initialstate.lower(3) = setup.bounds.phase(1).initialstate.lower(3) + auxdata.Stage2.mStruct*(mSPARTAN_vars(i)-1);
        setup_variations{i}.bounds.phase(1).initialstate.upper(3) = setup.bounds.phase(1).initialstate.upper(3) + auxdata.Stage2.mStruct*(mSPARTAN_vars(i)-1);
        setup_variations{i}.bounds.phase(1).state.lower(3) = setup.bounds.phase(1).state.lower(3) + auxdata.Stage2.mStruct*(mSPARTAN_vars(i)-1);
        setup_variations{i}.bounds.phase(1).state.upper(3) = setup.bounds.phase(1).state.upper(3) + auxdata.Stage2.mStruct*(mSPARTAN_vars(i)-1);
        setup_variations{i}.bounds.phase(1).finalstate.lower(3) = setup.bounds.phase(1).finalstate.lower(3) + auxdata.Stage2.mStruct*(mSPARTAN_vars(i)-1);
        setup_variations{i}.bounds.phase(1).finalstate.upper(3) = setup.bounds.phase(1).finalstate.upper(3) + auxdata.Stage2.mStruct*(mSPARTAN_vars(i)-1);
    
    
    end
elseif mode == 11 %
    for i = 1:length(mFuel_vars)  
        setup_variations{i} = setup;
        setup_variations{i}.bounds.phase(2).initialstate.lower(9) = Stage2.Initial.mFuel*mFuel_vars(i);
        setup_variations{i}.bounds.phase(2).initialstate.upper(9) = Stage2.Initial.mFuel*mFuel_vars(i);
        setup_variations{i}.bounds.phase(2).state.upper(9) = Stage2.Initial.mFuel*mFuel_vars(i);
        setup_variations{i}.bounds.phase(2).finalstate.upper(9) = Stage2.Initial.mFuel*mFuel_vars(i);
        
        setup_variations{i}.bounds.phase(1).initialstate.lower(3) = setup.bounds.phase(1).initialstate.lower(3) + Stage2.Initial.mFuel*(mFuel_vars(i)-1);
        setup_variations{i}.bounds.phase(1).initialstate.upper(3) = setup.bounds.phase(1).initialstate.upper(3) + Stage2.Initial.mFuel*(mFuel_vars(i)-1);
        setup_variations{i}.bounds.phase(1).state.lower(3) = setup.bounds.phase(1).state.lower(3) + Stage2.Initial.mFuel*(mFuel_vars(i)-1);
        setup_variations{i}.bounds.phase(1).state.upper(3) = setup.bounds.phase(1).state.upper(3) + Stage2.Initial.mFuel*(mFuel_vars(i)-1);
        setup_variations{i}.bounds.phase(1).finalstate.lower(3) = setup.bounds.phase(1).finalstate.lower(3) + Stage2.Initial.mFuel*(mFuel_vars(i)-1);
        setup_variations{i}.bounds.phase(1).finalstate.upper(3) = setup.bounds.phase(1).finalstate.upper(3) + Stage2.Initial.mFuel*(mFuel_vars(i)-1);
    
    end
end



for j = 1:length(setup_variations)
disp(['Starting Setup Variation ',num2str(j)])
    
for i = 1:num_it
setup_par(i) = setup_variations{j};
% setup_par(i).nlp.ipoptoptions.maxiterations = 500 + 10*i;
setup_par(i).guess.phase(2).state(:,1)   = [24000; 25000 + 1000*i]; % vary altitude guess
end

parfor i = 1:num_it
    
output_temp = gpops2(setup_par(i)); % Run GPOPS-2. Use setup for each parallel iteration.
% 
% % Run forward simulation for comparison between runs
% % Extract states
% alt22 = output_temp.result.solution.phase(3).state(:,1).';
% lon22 = output_temp.result.solution.phase(3).state(:,2).';
% lat22 = output_temp.result.solution.phase(3).state(:,3).';
% v22 = output_temp.result.solution.phase(3).state(:,4).'; 
% gamma22 = output_temp.result.solution.phase(3).state(:,5).'; 
% zeta22 = output_temp.result.solution.phase(3).state(:,6).';
% alpha22 = output_temp.result.solution.phase(3).state(:,7).';
% eta22 = output_temp.result.solution.phase(3).state(:,8).';
% mFuel22 = output_temp.result.solution.phase(3).state(:,9).'; 
% time22 = output_temp.result.solution.phase(3).time.';
% throttle22 = output_temp.result.solution.phase(3).state(:,10).';
% 
% 
% % Return Forward
% forward0 = [alt22(1),gamma22(1),v22(1),zeta22(1),lat22(1),lon22(1), mFuel22(1)];
% [f_t, f_y] = ode45(@(f_t,f_y) VehicleModelReturn_forward(f_t, f_y,auxdata,ControlInterp(time22,alpha22,f_t),ControlInterp(time22,eta22,f_t),ThrottleInterp(time22,throttle22,f_t)),time22,forward0);
% 
% error(i) = (f_y(end,6) + lon22(end))^2 + (f_y(end,5) + lat22(end))^2;
% error(i) = abs(mFuel22(end) - f_y(end,7));

input_test = output_temp.result.solution;
input_test.auxdata = auxdata;
phaseout_test = CombinedContinuous(input_test);
norm_error1 = [];
for num = [1:6 9]
% norm_error1(num,:) = (cumtrapz(output_temp.result.solution.phase(2).time,phaseout_test(2).dynamics(1:end,num)) + output_temp.result.solution.phase(2).state(1,num)- output_temp.result.solution.phase(2).state(:,num))./(max(output_temp.result.solution.phase(2).state(:,num))-min(output_temp.result.solution.phase(2).state(:,num)));
Stage2_int = cumtrapz([output_temp.result.solution.phase(2).time(1):0.1:output_temp.result.solution.phase(2).time(end)],pchip(output_temp.result.solution.phase(2).time,phaseout_test(2).dynamics(1:end,num),[output_temp.result.solution.phase(2).time(1):0.1:output_temp.result.solution.phase(2).time(end)]))';
norm_error1(num,:) = (interp1([output_temp.result.solution.phase(2).time(1):0.1:output_temp.result.solution.phase(2).time(end)],Stage2_int,output_temp.result.solution.phase(2).time)+ output_temp.result.solution.phase(2).state(1,num)- output_temp.result.solution.phase(2).state(:,num))./(max(output_temp.result.solution.phase(2).state(:,num))-min(output_temp.result.solution.phase(2).state(:,num)));

end
norm_error2 = [];
if returnMode == 1
for num = [1:6 9]
% norm_error2(num,:) = (cumtrapz(output_temp.result.solution.phase(3).time,phaseout_test(3).dynamics(1:end,num)) + output_temp.result.solution.phase(3).state(1,num)- output_temp.result.solution.phase(3).state(:,num))./(max(output_temp.result.solution.phase(3).state(:,num))-min(output_temp.result.solution.phase(3).state(:,num)));
Return_int = cumtrapz([output_temp.result.solution.phase(3).time(1):0.1:output_temp.result.solution.phase(3).time(end)],pchip(output_temp.result.solution.phase(3).time,phaseout_test(3).dynamics(1:end,num),[output_temp.result.solution.phase(3).time(1):0.1:output_temp.result.solution.phase(3).time(end)]))';
norm_error2(num,:) = (interp1([output_temp.result.solution.phase(3).time(1):0.1:output_temp.result.solution.phase(3).time(end)],Return_int,output_temp.result.solution.phase(3).time)+ output_temp.result.solution.phase(3).state(1,num)- output_temp.result.solution.phase(3).state(:,num))./(max(output_temp.result.solution.phase(3).state(:,num))-min(output_temp.result.solution.phase(3).state(:,num)));


end
else
 norm_error2 = 0;  
end
error(i) = max([max(abs(norm_error1)) max(abs(norm_error2))]);


PayloadMass(i) = -output_temp.result.objective;

output_store{i} = output_temp;

end

[min_error,index] = min(error); % Calculate the result which minimises the chosen error function

% [max_pl,index] = max(PayloadMass);% Calculate the result which maximises payload mass the chosen error function

output{j} = output_store{index};

% Clear output store to save memory and prevent write issues
if mode ~= 1
clear output_store
end

end

%% Process Results

Plotter(output,auxdata,auxdata.mode,auxdata.returnMode,auxdata.namelist,M_englist,T_englist,engine_data,MList_EngineOn,AOAList_EngineOn,mRocket,mSpartan,mFuel,h0,v0,bounds);









%% =========================================================================
% Troubleshooting Procedure
% =========================================================================
% 1: Check that you have posed your problem correctly ie. it is physically
% feasible and the bounds allow for a solution.
% 2: Check if there is any extrapolations causing bad dynamics. 
% 3: Check guess, is it reasonable?.
% 4: Increase number of iterations and collocation points.
% 5: Check for large nonlinearities, eg. atmospheric properties suddenly going to zero or thrust cutting off. These need to be eliminated or separated into phases.
% 6: Check for NaN values (check derivatives in Dynamics file while
% running).





function auxdata = AeroCalc(auxdata)


% Fetch aerodynamic data and compute interpolation splines.
% Calculate the flap deflection necessary for trim.
% Each set of aero corresponds to a different CG. 



Viscousaero_EngineOn = importdata('VC3D_viscousCoefficients_ascent.dat');
Viscousaero_EngineOff = importdata('VC3D_viscousCoefficients_descent.dat');

% These aerodynamic datasets have been created in ClicCalcCGVar.m

T_L = -1.327; % Thrust location in z, average (m), measured from CREO

% Full of fuel, with third stage
CG = 15.24;
CG_z = 0.049; % Centre of gravity location in z, calculated from Creo

aero1.aero_EngineOff = importdata('SPARTANaero15.24');
aero1.flapaero = importdata('SPARTANaeroFlaps15.24');
aero1.aero_EngineOn = importdata('SPARTANaeroEngineOn15.24');
aero1.aero_Engine = importdata('SPARTANEngine15.24');
aero1.Viscousaero_EngineOff = Viscousaero_EngineOff;
aero1.Viscousaero_EngineOn = Viscousaero_EngineOn;

[auxdata.interp.Cl_spline_EngineOff.fullFuel,auxdata.interp.Cd_spline_EngineOff.fullFuel,...
    auxdata.interp.Cl_spline_EngineOn.fullFuel,auxdata.interp.Cd_spline_EngineOn.fullFuel,...
    auxdata.interp.flap_spline_EngineOff.fullFuel,auxdata.interp.flap_spline_EngineOn.fullFuel,...
    auxdata.T_spline_Rear,auxdata.Fd_spline_NoEngine,auxdata.Cd_spline_ViscousEngineOff,auxdata.Cd_spline_ViscousEngineOn,...
    auxdata.L_spline_Rear,auxdata.T_spline] = AeroInt(aero1,auxdata,T_L,CG_z,CG);


% Cylindrical fuel tanks depleted, with third stage
CG = 15.13;
CG_z = 0.077;

aero2.aero_EngineOff = importdata('SPARTANaero15.13');
aero2.flapaero = importdata('SPARTANaeroFlaps15.13');
aero2.aero_EngineOn = importdata('SPARTANaeroEngineOn15.13');
aero2.aero_Engine = importdata('SPARTANEngine15.13');
aero2.Viscousaero_EngineOff = Viscousaero_EngineOff;
aero2.Viscousaero_EngineOn = Viscousaero_EngineOn;

[auxdata.interp.Cl_spline_EngineOff.cylTankEnd,auxdata.interp.Cd_spline_EngineOff.cylTankEnd,auxdata.interp.Cl_spline_EngineOn.cylTankEnd,auxdata.interp.Cd_spline_EngineOn.cylTankEnd,auxdata.interp.flap_spline_EngineOff.cylTankEnd,auxdata.interp.flap_spline_EngineOn.cylTankEnd] = AeroInt(aero2,auxdata,T_L,CG_z,CG);



% Empty, with third stage. 

CG = 15.74;
CG_z = 0.086;

aero3.aero_EngineOff = importdata('SPARTANaero15.74');
aero3.flapaero = importdata('SPARTANaeroFlaps15.74');
aero3.aero_EngineOn = importdata('SPARTANaeroEngineOn15.74');
aero3.aero_Engine = importdata('SPARTANEngine15.74');
aero3.Viscousaero_EngineOff = Viscousaero_EngineOff;
aero3.Viscousaero_EngineOn = Viscousaero_EngineOn;

[auxdata.interp.Cl_spline_EngineOff.noFuel,auxdata.interp.Cd_spline_EngineOff.noFuel,auxdata.interp.Cl_spline_EngineOn.noFuel,auxdata.interp.Cd_spline_EngineOn.noFuel,auxdata.interp.flap_spline_EngineOff.noFuel,auxdata.interp.flap_spline_EngineOn.noFuel] = AeroInt(aero3,auxdata,T_L,CG_z,CG);

% Flyback, empty without third stage, empty.
CG = 15.16;
CG_z = -0.2134; % calculated fom CREO

aero4.aero_EngineOff = importdata('SPARTANaero15.16');
aero4.flapaero = importdata('SPARTANaeroFlaps15.16');
aero4.aero_EngineOn = importdata('SPARTANaeroEngineOn15.16');
aero4.aero_Engine = importdata('SPARTANEngine15.16');
aero4.Viscousaero_EngineOff = Viscousaero_EngineOff;
aero4.Viscousaero_EngineOn = Viscousaero_EngineOn;

[auxdata.interp.Cl_spline_EngineOff.noThirdStageEmpty,...
    auxdata.interp.Cd_spline_EngineOff.noThirdStageEmpty,auxdata.interp.Cl_spline_EngineOn.noThirdStageEmpty,...
    auxdata.interp.Cd_spline_EngineOn.noThirdStageEmpty,auxdata.interp.flap_spline_EngineOff.noThirdStageEmpty,...
    auxdata.interp.flap_spline_EngineOn.noThirdStageEmpty] = AeroInt(aero4,auxdata,T_L,CG_z,CG);

% Flyback, empty without third stage, conical fuel tank full.
CG = 14.3;
CG_z = -0.183; % calculated fom CREO

aero5.aero_EngineOff = importdata('SPARTANaero14.3');
aero5.flapaero = importdata('SPARTANaeroFlaps14.3');
aero5.aero_EngineOn = importdata('SPARTANaeroEngineOn14.3');
aero5.aero_Engine = importdata('SPARTANEngine14.3');
aero5.Viscousaero_EngineOff = Viscousaero_EngineOff;
aero5.Viscousaero_EngineOn = Viscousaero_EngineOn;

[auxdata.interp.Cl_spline_EngineOff.noThirdStagecylTankEnd,...
    auxdata.interp.Cd_spline_EngineOff.noThirdStagecylTankEnd,auxdata.interp.Cl_spline_EngineOn.noThirdStagecylTankEnd,...
    auxdata.interp.Cd_spline_EngineOn.noThirdStagecylTankEnd,auxdata.interp.flap_spline_EngineOff.noThirdStagecylTankEnd,...
    auxdata.interp.flap_spline_EngineOn.noThirdStagecylTankEnd] = AeroInt(aero5,auxdata,T_L,CG_z,CG);


clear aero1
clear aero2
clear aero3
clear aero4
clear aero5
clear Viscousaero_EngineOn
clear Viscousaero_EngineOff
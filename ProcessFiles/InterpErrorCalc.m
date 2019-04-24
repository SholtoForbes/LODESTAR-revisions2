%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to estimate the error in the interpolation of aerodynamic and
% atmospheric conditions.
clear all;
clc



%% Add Necessary Paths


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('..\EngineData')
addpath('..\')
addpath('.\Interpolators')
addpath('.\Processing')
addpath('.\Forward Simulations')
addpath('.\Dynamics')
addpath ('..\ViscousAero')
addpath('..\CG15.16') % Add folder containing aerodynamic datasets


%% FIRST STAGE


%  auxdata.Throttle = .7 % throttle the Merlin engine down by a modeant value, to enable easier pitchover

% Aerodynamics File Path
Aero1_Full = dlmread('FirstStageAero23.365');
Aero1_Empty = dlmread('FirstStageAero16.8631');

%% Vehicle 

% mRocket =21816; % total mass of scaled Falcon at 9.5m, note, this will not be the final total mass. Calculated using the method outlined in SIZING.docx
% mRocket =19569; % total mass of scaled Falcon at 8.5m, note, this will
% not be the final total mass. Calculated using the method outlined in
% SIZING.docx % these use wrong merlin mass
% mRocket =18536; % total mass of scaled Falcon at 8m, leaving 0.5m for the merlin engine, note, this may not be the final total mass. Calculated using the method outlined in SIZING.docx
% mRocket =17417; % total mass of scaled Falcon at 7.5m, leaving 1m for the merlin engine, note, this may not be the final total mass. Calculated using the method outlined in SIZING.docx
mRocket =18088.5; % total mass of scaled Falcon at 7.8m, leaving 0.7m for the merlin engine (measured in users guide from tip to gimbal end), note, this may not be the final total mass. Calculated using the method outlined in SIZING.docx

auxdata.Vehicle.mRocket = mRocket;
mEngine = 630; % Mass of Merlin 1C
mFuel = 0.953*(mRocket-mEngine); % structural mass fraction calculated without engine
auxdata.Vehicle.mFuel=mFuel;
mSpartan = 9819.11;

% Thrust and Isp are modified with altitude through the formula: %
% SL + (101325-P_atm)*Mod %
 
auxdata.Vehicle.T.SL = 555900; % Thrust from Falcon 1 users guide. 
auxdata.Vehicle.T.Mod = 0.5518; % exit area calculated in SIZING.docx

auxdata.Vehicle.Isp.SL = 275; % linear regression of SL and vacuum Isp. From encyclopaedia astronautica, backed up by falcon 1 users guide
auxdata.Vehicle.Isp.Mod = 2.9410e-04;

auxdata.Vehicle.Area = 62.77; 

auxdata.Vehicle.CG_Full = 23.8257; % centre of gravity (m) calculated in clicalcCGvar.m
auxdata.Vehicle.CG_Empty = 16.7507; % centre of gravity (m) calculated in clicalcCGvar.m

auxdata.Vehicle.L = 22.94+8.5;
%% Import Atmosphere

auxdata.Atmosphere = dlmread('atmosphere.txt');

%% Calculate Aerodynamic Splines

interp.Lift = scatteredInterpolant(Aero1_Full(:,1),Aero1_Full(:,2),Aero1_Full(:,3));
interp.Drag = scatteredInterpolant(Aero1_Full(:,1),Aero1_Full(:,2),Aero1_Full(:,4));
interp.MomentFull = scatteredInterpolant(Aero1_Full(:,1),Aero1_Full(:,2),Aero1_Full(:,5));
interp.MomentEmpty = scatteredInterpolant(Aero1_Empty(:,1),Aero1_Empty(:,2),Aero1_Empty(:,5));

M_list = unique(sort(Aero1_Full(:,1))); % create unique list of Mach numbers from engine data
M_interp = unique(sort(Aero1_Full(:,1)));

AoA_list = unique(sort(Aero1_Full(:,2))); % create unique list of angle of attack numbers from engine data 
AoA_interp = unique(sort(Aero1_Full(:,2)));

[grid.M,grid.AoA] =  ndgrid(M_interp,AoA_interp);
grid.Lift = interp.Lift(grid.M,grid.AoA);
grid.Drag = interp.Drag(grid.M,grid.AoA);
grid.MomentFull = interp.MomentFull(grid.M,grid.AoA);
grid.MomentEmpty = interp.MomentEmpty(grid.M,grid.AoA);

auxdata.interp.LiftGridded = griddedInterpolant(grid.M,grid.AoA,grid.Lift,'spline','linear');
auxdata.interp.DragGridded = griddedInterpolant(grid.M,grid.AoA,grid.Drag,'spline','linear');
auxdata.interp.MomentGriddedFull = griddedInterpolant(grid.M,grid.AoA,grid.MomentFull,'spline','linear');
auxdata.interp.MomentGriddedEmpty = griddedInterpolant(grid.M,grid.AoA,grid.MomentEmpty,'spline','linear');

%% Create Viscous Drag Interpolant

Aero_Visc1 = dlmread('VC3D_C_v_SPARTANFALCON_descent.dat');


interp.Cl_scattered_Viscousaero_1 = scatteredInterpolant(Aero_Visc1(:,1),Aero_Visc1(:,2),Aero_Visc1(:,3),Aero_Visc1(:,4));
interp.Cd_scattered_Viscousaero_1 = scatteredInterpolant(Aero_Visc1(:,1),Aero_Visc1(:,2),Aero_Visc1(:,3),Aero_Visc1(:,5));

% MList_Visc1 = [0; unique(Aero_Visc1(:,1))]; % 0 included to prevent bad extrapolation
% 
% 
% AoAList_Visc1 = unique(Aero_Visc1(:,2));
% 
% altList_Visc1 = unique(Aero_Visc1(:,3)); 

MList_Visc1 = [0:0.1: max(unique(Aero_Visc1(:,1)))]; % 0 included to prevent bad extrapolation


AoAList_Visc1 = [min(unique(Aero_Visc1(:,2))):0.1:max(unique(Aero_Visc1(:,2)))];

altList_Visc1 = [min(unique(Aero_Visc1(:,3))):500:max(unique(Aero_Visc1(:,3)))];

[Mgrid_Visc1,AOAgrid_Visc1,altgrid_Visc1] = ndgrid(MList_Visc1,AoAList_Visc1,altList_Visc1);

Cl_Grid_Visc1 = [];
Cd_Grid_Visc1 = [];

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

auxdata.interp.Cl_gridded_Visc1 = griddedInterpolant(Mgrid_Visc1,AOAgrid_Visc1,altgrid_Visc1,Cl_Grid_Visc1,'spline','linear'); % Note: spline interp here was found to produce some odd results in in-between values if mesh of points was not refined as above
auxdata.interp.Cd_gridded_Visc1 = griddedInterpolant(Mgrid_Visc1,AOAgrid_Visc1,altgrid_Visc1,Cd_Grid_Visc1,'spline','linear');




%% Atmosphere Data %%======================================================
% Fetch atmospheric data and compute interpolation splines.

Atmosphere = dlmread('atmosphere.txt');
interp.Atmosphere = Atmosphere;
auxdata.interp.Atmosphere = interp.Atmosphere;

auxdata.interp.c_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,5)); % Calculate speed of sound using atmospheric data
auxdata.interp.rho_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,4)); % Calculate density using atmospheric data
auxdata.interp.mu_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,6)); % Calculate dynamic viscosity using atmospheric data
auxdata.interp.p_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,3)); % Calculate pressure using atmospheric data
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
auxdata.A_flap = 2.046;
auxdata.l_flap = 1.55;
auxdata.Flap_Base = 19.28; %Distance from nose tip to base (connected to wings) of flaps
auxdata.L = 22.94;
%% Third Stage Aerodynamic Data

% auxdata.Aero3 = dlmread('ThirdStageAeroCoeffs.txt');
% 
% Aero3 = auxdata.Aero3;
% auxdata.interp.Drag_interp3 = scatteredInterpolant(Aero3(:,1),Aero3(:,2),Aero3(:,5));
% % 
% auxdata.interp.Lift_interp3 = scatteredInterpolant(Aero3(:,1),Aero3(:,2),Aero3(:,6));
% 
% auxdata.interp.CP_interp3 = scatteredInterpolant(Aero3(:,1),Aero3(:,2),Aero3(:,7));
% 
% auxdata.interp.CN_interp3 = scatteredInterpolant(Aero3(:,1),Aero3(:,2),Aero3(:,4));
% 
% auxdata.interp.Max_AoA_interp3 = scatteredInterpolant(Aero3(:,1),Aero3(:,4),Aero3(:,2));

aero30km = datcomimport('for006-30km.dat'); % read missile datcom data files
aero50km = datcomimport('for006-50km.dat');
aero70km = datcomimport('for006-70km.dat');
aero90km = datcomimport('for006-90km.dat');

altlist3 = [30000 50000 70000 90000]; 
Mlist3 = aero30km{1}.mach;
AoAlist3 = aero30km{1}.alpha;

%Drag
[Mgrid3,AoAgrid3,altgrid3] = ndgrid(Mlist3,AoAlist3,altlist3);

cdgrid3 = aero30km{1}.cd;
cdgrid3(:,:,2) = aero50km{1}.cd;
cdgrid3(:,:,3) = aero70km{1}.cd;
cdgrid3(:,:,4) = aero90km{1}.cd;

P = [2 1 3]; %permute to ndgrid format
cdgrid3 = permute(cdgrid3, P);

auxdata.interp.Drag_interp3 = griddedInterpolant(Mgrid3,AoAgrid3,altgrid3,cdgrid3,'spline','linear');

%Lift
[Mgrid3,AoAgrid3,altgrid3] = ndgrid(Mlist3,AoAlist3,altlist3);

clgrid3 = aero30km{1}.cl;
clgrid3(:,:,2) = aero50km{1}.cl;
clgrid3(:,:,3) = aero70km{1}.cl;
clgrid3(:,:,4) = aero90km{1}.cl;

P = [2 1 3]; %permute to ndgrid format
clgrid3 = permute(clgrid3, P);

auxdata.interp.Lift_interp3 = griddedInterpolant(Mgrid3,AoAgrid3,altgrid3,clgrid3,'spline','linear');

%Normal Force
[Mgrid3,AoAgrid3,altgrid3] = ndgrid(Mlist3,AoAlist3,altlist3);

cngrid3 = aero30km{1}.cn;
cngrid3(:,:,2) = aero50km{1}.cn;
cngrid3(:,:,3) = aero70km{1}.cn;
cngrid3(:,:,4) = aero90km{1}.cn;

P = [2 1 3]; %permute to ndgrid format
cngrid3 = permute(cngrid3, P);

auxdata.interp.CN_interp3 = griddedInterpolant(Mgrid3,AoAgrid3,altgrid3,cngrid3,'spline','linear');

%CP Location
[Mgrid3,AoAgrid3,altgrid3] = ndgrid(Mlist3,AoAlist3,altlist3);

xcpgrid3 = aero30km{1}.xcp;
xcpgrid3(:,:,2) = aero50km{1}.xcp;
xcpgrid3(:,:,3) = aero70km{1}.xcp;
xcpgrid3(:,:,4) = aero90km{1}.xcp;

P = [2 1 3]; %permute to ndgrid format
xcpgrid3 = permute(xcpgrid3, P);

auxdata.interp.CP_interp3 = griddedInterpolant(Mgrid3,AoAgrid3,altgrid3,xcpgrid3,'spline','linear');


%Max Aoa (reverse interp with CN)
Aero3forMaxAoa = dlmread('ThirdStageAeroCoeffs.txt'); % import aerodynamics with approximated altitudes for max aoa calculations
auxdata.interp.Max_AoA_interp3 = scatteredInterpolant(Aero3forMaxAoa (:,1),Aero3forMaxAoa (:,4),Aero3forMaxAoa (:,2));

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


%%
% Created by Sholto Forbes 25/7/18
% This file processes the GPOPS-2 solution output for the
% rocket-scramjet-rocket launch trajectory.
%%

function [] = Plotter(output,auxdata,mode,returnMode,namelist,M_englist,T_englist,engine_data,MList_EngineOn,AOAList_EngineOn,mRocket,mSpartan,mFuel,h0,v0,bounds)

addpath('.\Processing\addaxis')
addpath('.\Processing\axlabel')
addpath('.\Processing\num2words')

Timestamp = datestr(now,30)
mkdir('../ArchivedResults', strcat(Timestamp, 'mode', num2str(mode),num2str(returnMode)))


copyfile('CombinedProbGPOPS.m',sprintf('../ArchivedResults/%s/SecondStageProb.m',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))))

save output
movefile('output.mat',sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));


% movefile(strcat('SPARTAN-CombinedIPOPTinfo.txt'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));


% % Run multiple times if output has mutltiple cells
if mode ~= 99
for j = 1:length(output)


auxdata = output{j}.result.setup.auxdata

% =========================================================================
% Assign the primal variables

time1 = output{j}.result.solution.phase(1).time.';
% 
alt1 = output{j}.result.solution.phase(1).state(:,1).';
v1 = output{j}.result.solution.phase(1).state(:,2).';
m1 = output{j}.result.solution.phase(1).state(:,3).';
gamma1 = output{j}.result.solution.phase(1).state(:,4).';
alpha1 = output{j}.result.solution.phase(1).state(:,5).';
zeta1 = output{j}.result.solution.phase(1).state(:,6).';
lat1 = output{j}.result.solution.phase(1).state(:,8).';
lon1 = output{j}.result.solution.phase(1).state(:,9).';

Throttle1 = output{j}.result.solution.phase(1).state(:,10).';
% 

alt21 = output{j}.result.solution.phase(2).state(:,1);
alt22 = output{j}.result.solution.phase(3).state(:,1);
lon21 = output{j}.result.solution.phase(2).state(:,2);
lon22 = output{j}.result.solution.phase(3).state(:,2);
lat21 = output{j}.result.solution.phase(2).state(:,3);
lat22 = output{j}.result.solution.phase(3).state(:,3);
v21 = output{j}.result.solution.phase(2).state(:,4); 
v22 = output{j}.result.solution.phase(3).state(:,4); 
gamma21 = output{j}.result.solution.phase(2).state(:,5); 
gamma22 = output{j}.result.solution.phase(3).state(:,5); 
zeta21 = output{j}.result.solution.phase(2).state(:,6);
zeta22 = output{j}.result.solution.phase(3).state(:,6);
alpha21 = output{j}.result.solution.phase(2).state(:,7);
alpha22 = output{j}.result.solution.phase(3).state(:,7);
eta21 = output{j}.result.solution.phase(2).state(:,8);
eta22 = output{j}.result.solution.phase(3).state(:,8);
mFuel21 = output{j}.result.solution.phase(2).state(:,9); 
mFuel22 = output{j}.result.solution.phase(3).state(:,9); 

throttle22 = output{j}.result.solution.phase(3).state(:,10);

aoadot21  = output{j}.result.solution.phase(2).control(:,1); 
etadot21  = output{j}.result.solution.phase(2).control(:,2); 

aoadot22  = output{j}.result.solution.phase(3).control(:,1); 
etadot22  = output{j}.result.solution.phase(3).control(:,2); 

time21 = output{j}.result.solution.phase(2).time;
time22 = output{j}.result.solution.phase(3).time;


alt3  = output{j}.result.solution.phase(4).state(:,1);
v3    = output{j}.result.solution.phase(4).state(:,2);
gamma3  = output{j}.result.solution.phase(4).state(:,3);
m3    = output{j}.result.solution.phase(4).state(:,4);
aoa3    = output{j}.result.solution.phase(4).state(:,5);
lat3    = output{j}.result.solution.phase(4).state(:,6);
zeta3    = output{j}.result.solution.phase(4).state(:,7);
aoadot3       = output{j}.result.solution.phase(4).control(:,1);

time3 = output{j}.result.solution.phase(4).time;


figure(01)
subplot(9,1,1)
hold on
plot(time21,alt21)
plot(time22,alt22)
subplot(9,1,2)
hold on
plot(time21,v21)
plot(time22,v22)
subplot(9,1,3)
hold on
plot(time21,lon21)
plot(time22,lon22)
subplot(9,1,4)
hold on
plot(time21,lat21)
plot(time22,lat22)
subplot(9,1,5)
hold on
plot(time21,v21)
plot(time22,v22)
subplot(9,1,6)
hold on
plot(time21,gamma21)
plot(time22,gamma22)
subplot(9,1,7)
hold on
plot(time21,ones(1,length(time21)))
plot(time22,throttle22)


figure(230)
hold on
plot3(lon21,lat21,alt21)
plot3(lon22,lat22,alt22)

  
%     figure(2301)
% hold on
% 
% axesm('pcarree','Origin',[0 rad2deg(lon21(1)) 0])
% geoshow('landareas.shp','FaceColor',[0.8 .8 0.8])
% % plotm(rad2deg(lat),rad2deg(lon+lon0))
% plotm(rad2deg(lat21),rad2deg(lon21),'b')
% plotm(rad2deg(lat22),rad2deg(lon22),'r')
%     
%     cities = shaperead('worldcities', 'UseGeoCoords', true);
% lats = extractfield(cities,'Lat');
% lons = extractfield(cities,'Lon');
% geoshow(lats, lons,...
%         'DisplayType', 'point',...
%         'Marker', 'o',...
%         'MarkerEdgeColor', 'r',...
%         'MarkerFaceColor', 'r',...
%         'MarkerSize', 2)

% =========================================================================

%% Third Stage
% simulate unpowered third stage trajectory from end point

[rdot3,xidot3,phidot3,gammadot3,vdot3,zetadot3, mdot3, Vec_angle3, AoA_max3, T3, L3, D3, q3] = ThirdStageDynamics(alt3,gamma3,v3,m3,aoa3,time3,auxdata,aoadot3,lat3,zeta3);

figure(3223)
plot(time3,rad2deg(aoa3))
hold on
plot(time3,rad2deg(AoA_max3))
xlabel('Time (s)')
ylabel('Angle (deg)')
legend('AoA','AoA_max')
ylim([0 20])

lon3 = [];
lon3(1) = lon21(end);
for i = 2:length(time3)
    lon3(i) = lon3(i-1) + xidot3(i-1)*(time3(i)-time3(i-1));
end

[AltF_actual, v3F, altexo, v3exo, timeexo, mpayload, Alpha3, mexo,qexo,gammaexo,Dexo,zetaexo,latexo,incexo,Texo,CLexo,Lexo,incdiffexo,lonexo,dvtot3,m3_4,v3exo_coordchange,hs] = ThirdStageSim(alt3(end),gamma3(end),v3(end), lat3(end),lon3(end), zeta3(end), m3(end), auxdata);

[~, ~, ~, ~, ~, mpayloadtest] = ThirdStageSim(alt3(end),gamma3(end),v3(end), lat3(end)-0.008,lon3(end), zeta3(end)+0.0005, m3(end), auxdata);

payload_latdiff = mpayload-mpayloadtest

% ThirdStagePayloadMass = -output{j}.result.objective;

ThirdStagePayloadMass = mpayload;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          OUTPUT             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nodes = length(alt21)

[altdot21,xidot21,phidot21,gammadot21,a21,zetadot21, q21, M21, Fd21, rho021,L21,Fueldt21,T21,Isp21,q121,flapdeflection21,heating_rate21,CG21,T121,P121,M121,P021,T021,P_1_tip21,T_1_tip21,rho_1_tip21,M_1_tip21,v_1_tip21,heating_rate_LE21] = SPARTANDynamics(gamma21, alt21, v21,auxdata,zeta21,lat21,lon21,alpha21,eta21,1.*ones(length(alt21),1), mFuel21,mFuel21(1),mFuel21(end), 1, 0);

[~,~,~,~,~,~, q22, M22, Fd22, rho022,L22,Fueldt22,T22,Isp22,q122,flapdeflection22,heating_rate22,CG22,T122,P122,M122,P022,T022,P_1_tip22,T_1_tip22,rho_1_tip22,M_1_tip22,v_1_tip22,heating_rate_LE22] = SPARTANDynamics(gamma22, alt22, v22,auxdata,zeta22,lat22,lon22,alpha22,eta22,throttle22, mFuel22,0,0, 0, 0);

%% Modify throttle for return 


throttle22(M22<5.0) = throttle22(M22<5.0).*gaussmf(M22(M22<5.0),[.01,5]); % remove nonsense throttle points
throttle22(q122<20000) = throttle22(q122<20000).*gaussmf(q122(q122<20000),[100,20000]); % rapidly reduce throttle to 0 after passing the lower limit of 20kPa dynamic pressure. This dynamic pressure is after the conical shock.
    
Isp22(M22<5.0) = Isp22(M22<5.0).*gaussmf(M22(M22<5.0),[.01,5]); %
Isp22(q122<20000) = Isp22(q122<20000).*gaussmf(M22(q122<20000),[.1,5]); %


% Separation_LD = lift(end)/Fd(end)

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% FORWARD SIMULATION
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% This is a full forward simulation, using the angle of attack and flap
% deflection at each node.

% Note, because the nodes are spaced widely, small interpolation
% differences result in the forward simulation being slightly different
% than the actual. This is mostly a check to see if they are close. 


forward0 = [alt21(1),gamma21(1),v21(1),zeta21(1),lat21(1),lon21(1), mFuel21(1)];

[f_t, f_y] = ode45(@(f_t,f_y) VehicleModelAscent_forward(f_t, f_y,auxdata,ControlInterp(time21,alpha21,f_t),ControlInterp(time21,eta21,f_t),1,mFuel21(1),mFuel21(end)),time21(1:end),forward0);

forward_error21 = [(alt21-f_y(:,1))/(max(alt21)-min(alt21)) (lon21-f_y(:,6))/(max(lon21)-min(lon21))...
    (lat21-f_y(:,5))/(max(lat21)-min(lat21)) (v21-f_y(:,3))/(max(v21)-min(v21))...
    (gamma21-f_y(:,2))/(max(gamma21)-min(gamma21)) (zeta21-f_y(:,4))/(max(zeta21)-min(zeta21)) (mFuel21-f_y(:,7))/(max(mFuel21)-min(mFuel21))];



figure(212)
subplot(7,1,[1 2])
hold on
plot(f_t(1:end),f_y(:,1));
plot(time21,alt21);

subplot(7,1,3)
hold on
plot(f_t(1:end),f_y(:,2));
plot(time21,gamma21);


subplot(7,1,4)
hold on
plot(f_t(1:end),f_y(:,3));
plot(time21,v21);

subplot(7,1,6)
hold on
plot(f_t(1:end),f_y(:,4));
plot(time21,zeta21);

subplot(7,1,7)
hold on
plot(f_t(1:end),f_y(:,7));
plot(time21,mFuel21);



% Return Forward
forward0 = [alt22(1),gamma22(1),v22(1),zeta22(1),lat22(1),lon22(1), mFuel22(1)];
% 
% [f_t, f_y] = ode45(@(f_t,f_y) VehicleModelReturn_forward(f_t, f_y,auxdata,ControlInterp(time22,alpha22,f_t),ControlInterp(time22,eta22,f_t),ThrottleInterp(time22,throttle22,f_t)),time22(1):time22(end),forward0);
% [f_t, f_y] = ode45(@(f_t,f_y) VehicleModelReturn_forward(f_t, f_y,auxdata,ControlInterp(time22,alpha22,f_t),ControlInterp(time22,eta22,f_t),ThrottleInterp(time22,throttle22,f_t)),time22,forward0);


[f_t221, f_y221] = ode15s(@(f_t,f_y) VehicleModelReturn_forward(f_t, f_y,auxdata,ControlInterp(time22,alpha22,f_t),ControlInterp(time22,eta22,f_t),ThrottleInterp(time22,throttle22,f_t)),time22(1:ceil(length(time22)/6)-1),forward0);

forward0 = [alt22(ceil(length(time22)/6)),gamma22(ceil(length(time22)/6)),v22(ceil(length(time22)/6)),zeta22(ceil(length(time22)/6)),lat22(ceil(length(time22)/6)),lon22(ceil(length(time22)/6)), mFuel22(ceil(length(time22)/6))];

[f_t222, f_y222] = ode15s(@(f_t,f_y) VehicleModelReturn_forward(f_t, f_y,auxdata,ControlInterp(time22,alpha22,f_t),ControlInterp(time22,eta22,f_t),ThrottleInterp(time22,throttle22,f_t)),time22(ceil(length(time22)/6):ceil(length(time22)/3)-1),forward0);


forward0 = [alt22(ceil(length(time22)/3)),gamma22(ceil(length(time22)/3)),v22(ceil(length(time22)/3)),zeta22(ceil(length(time22)/3)),lat22(ceil(length(time22)/3)),lon22(ceil(length(time22)/3)), mFuel22(ceil(length(time22)/3))];

[f_t223, f_y223] = ode15s(@(f_t,f_y) VehicleModelReturn_forward(f_t, f_y,auxdata,ControlInterp(time22,alpha22,f_t),ControlInterp(time22,eta22,f_t),ThrottleInterp(time22,throttle22,f_t)),time22(ceil(length(time22)/3):end),forward0);


f_t2 = [f_t221 ; f_t222 ; f_t223];
% 
f_y2 = [f_y221; f_y222 ; f_y223];

forward_error22 = [(alt22-f_y2(:,1))/(max(alt22)-min(alt22)) (lon22-f_y2(:,6))/(max(lon22)-min(lon22))...
    (lat22-f_y2(:,5))/(max(lat22)-min(lat22)) (v22-f_y2(:,3))/(max(v22)-min(v22))...
    (gamma22-f_y2(:,2))/(max(gamma22)-min(gamma22)) (zeta22-f_y2(:,4))/(max(zeta22)-min(zeta22)) (mFuel22-f_y2(:,7))/(max(mFuel22)-min(mFuel22))];


figure(213)
subplot(7,1,1)
hold on
plot(f_t2(1:end),f_y2(:,1));
plot(time22,alt22);

% gamma  = output.result.solution.phase.state(:,5);

subplot(7,1,2)
hold on
plot(f_t2(1:end),f_y2(:,2));
plot(time22,gamma22);

% latitude  = output.result.solution.phase.state(:,3);
subplot(7,1,3:5)
hold on
plot(f_y2(:,6),f_y2(:,5));
plot(lon22,lat22);

subplot(7,1,6)
hold on
plot(f_t2(1:end),f_y2(:,7));
plot(time22,mFuel22);



%% plot Ascent
addpath('addaxis')
addpath('axlabel')

 figure(211)
 fig = gcf;
set(fig,'Position',[200 0 850 1200])
dim = [.62 .78 .2 .2];
annotation('textbox',dim,'string',{['Payload Mass: ', num2str(ThirdStagePayloadMass), ' kg'],['Fuel Used: ' num2str(1562 - mFuel21(end)) ' kg']},'FitBoxToText','on');  
hold on
% %suptitle('Second Stage Ascent')
subplot(4,1,1)
% set(gca,'xticklabels',[])
hold on
xlim([0 time21(end)-time21(1)]);
 plot(time21-time21(1),alt21/1000,'-','color','k', 'linewidth', 1.);
% ylim([-30,40]);
ylabel('altitude(km)');
addaxis(time21-time21(1),M21,'--','color','k', 'linewidth', 1.);
addaxislabel(2,'Mach no.');

% addaxis(time,fpa,':','color','k', 'linewidth', 1.);
% addaxislabel(3,'Trajectory Angle (deg)');


addaxis(time21-time21(1),q21/1000,':','color','k', 'linewidth', 1.2);
addaxislabel(3,'Dynamic Pressure (kPa)');

legend(  'Altitude', 'Mach no.', 'Dynamic Pressure', 'location', 'best');


subplot(4,1,2)
% set(gca,'xticklabels',[])
hold on
xlim([0 time21(end)-time21(1)]);
 plot(time21-time21(1),rad2deg(gamma21),'-','color','k', 'linewidth', 1.);
% ylim([-30,40]);
ylabel('Trajectory Angle (deg)');

addaxis(time21-time21(1),rad2deg(zeta21),'--','color','k', 'linewidth', 1.);
addaxislabel(2,'Heading Angle (deg)');


LD21 = (L21 - auxdata.L_spline_Rear(M21,rad2deg(alpha21),alt21/1000)*auxdata.A.*q21)./(Fd21 - auxdata.T_spline_Rear(M21,rad2deg(alpha21),alt21/1000)*auxdata.A.*q21);

addaxis(time21-time21(1),LD21,':','color','k', 'linewidth', 1.2);
addaxislabel(3,'L/D');

legend(  'Trajectory Angle', 'Heading Angle', 'L/D', 'location', 'best');


subplot(4,1,3)
% set(gca,'xticklabels',[])
hold on
xlim([0 time21(end)-time21(1)]);
 plot(time21-time21(1),rad2deg(alpha21),'-','color','k', 'linewidth', 1.);
% ylim([-30,40]);
ylabel('Angle of Attack (deg)');

if returnMode == 1
addaxis(time21-time21(1),rad2deg(eta21),'--','color','k', 'linewidth', 1.);
addaxislabel(2,'Bank Angle (deg)');
else
   addaxis(time21-time21(1),CG21,'--','color','k', 'linewidth', 0.1);
addaxislabel(2,'Centre of Gravity'); 
end

addaxis(time21-time21(1),flapdeflection21,':','color','k', 'linewidth', 1.2);
addaxislabel(3,'Flap Deflection (deg)');

if returnMode == 1
legend(  'Angle of Attack', 'Bank Angle', 'Flap Deflection', 'location', 'best');
else
 legend(  'Angle of Attack', 'Centre of Gravity', 'Flap Deflection', 'location', 'best');
end

subplot(4,1,4)
Isp21 = T21./Fueldt21./9.81;
IspNet21 = (T21-Fd21)./Fueldt21./9.81;
xlim([0 time21(end)-time21(1)]);
hold on

 plot(time21-time21(1),IspNet21,'-','color','k', 'linewidth', 1.);
% ylim([-30,40]);
ylabel('Net Isp (s)');
xlabel('Time (s)');

addaxis(time21-time21(1),T21/1000,'--','color','k', 'linewidth', 1.);
addaxislabel(2,'Thrust (kN)');

addaxis(time21-time21(1),mFuel21,':','color','k', 'linewidth', 1.);
addaxislabel(3,'Fuel Mass (kg)');

legend(  'Net Isp', 'Thrust', 'Fuel Mass', 'location', 'best');


figure(221)
 fig = gcf;
set(fig,'Position',[200 0 850 1200])
dim = [.62 .78 .2 .2];
annotation('textbox',dim,'string',{['Fuel Used: ' num2str(mFuel22(1)) ' kg']},'FitBoxToText','on');  
hold on
% %suptitle('Second Stage Ascent')
subplot(4,1,1)
% set(gca,'xticklabels',[])
hold on
xlim([0 time22(end)-time22(1)]);
 plot(time22-time22(1),alt22/1000,'-','color','k', 'linewidth', 1.);
% ylim([-30,40]);
ylabel('altitude(km)');
addaxis(time22-time22(1),M22,'--','color','k', 'linewidth', 1.);
addaxislabel(2,'Mach no.');

% addaxis(time,fpa,':','color','k', 'linewidth', 1.);
% addaxislabel(3,'Trajectory Angle (deg)');


addaxis(time22-time22(1),q22/1000,':','color','k', 'linewidth', 1.2);
addaxislabel(3,'Dynamic Pressure (kPa)');

legend(  'Altitude', 'Mach no.', 'Dynamic Pressure', 'location', 'best');


subplot(4,1,2)
% set(gca,'xticklabels',[])
hold on
xlim([0 time22(end)-time22(1)]);
 plot(time22-time22(1),rad2deg(gamma22),'-','color','k', 'linewidth', 1.);
% ylim([-30,40]);
ylabel('Trajectory Angle (deg)');

addaxis(time22-time22(1),rad2deg(zeta22),'--','color','k', 'linewidth', 1.);
addaxislabel(2,'Heading Angle (deg)');

LD22 = (L22 - auxdata.L_spline_Rear(M22,rad2deg(alpha22),alt22/1000)*auxdata.A.*q22.*throttle22)./(Fd22 - auxdata.T_spline_Rear(M22,rad2deg(alpha22),alt22/1000)*auxdata.A.*q22.*throttle22);

addaxis(time22-time22(1),LD22,':','color','k', 'linewidth', 1.2);
addaxislabel(3,'L/D');

legend(  'Trajectory Angle', 'Heading Angle', 'L/D', 'location', 'best');


subplot(4,1,3)
% set(gca,'xticklabels',[])
hold on
xlim([0 time22(end)-time22(1)]);
 plot(time22-time22(1),rad2deg(alpha22),'-','color','k', 'linewidth', 1.);
% ylim([-30,40]);
ylabel('Angle of Attack (deg)');

addaxis(time22-time22(1),rad2deg(eta22),'--','color','k', 'linewidth', 1.);
addaxislabel(2,'Bank Angle (deg)');

addaxis(time22-time22(1),flapdeflection22,':','color','k', 'linewidth', 1.2);
addaxislabel(3,'Flap Deflection (deg)');

legend(  'Angle of Attack', 'Bank Angle', 'Flap Deflection', 'location', 'best');

subplot(4,1,4)

xlim([0 time22(end)-time22(1)]);
hold on

 plot(time22-time22(1),Isp22,'-','color','k', 'linewidth', 1.);
% ylim([-30,40]);
ylabel('Potential Isp (s)');
xlabel('Time (s)');

addaxis(time22-time22(1),T22/1000,'--','color','k', 'linewidth', 1.);
addaxislabel(2,'Thrust (kN)');

addaxis(time22-time22(1),throttle22,':','color','k', 'linewidth', 1.);
addaxislabel(3,'Throttle');

legend(  'Net Isp', 'Thrust', 'Throttle', 'location', 'best');



ThirdStageFlag = [ones(length(time21),1); zeros(length(time22),1)];

SecondStageStates = [[time21; time22] [alt21; alt22] [M21;M22] [lon21; lon22] [lat21; lat22] [v21; v22] [gamma21; gamma22] [zeta21; zeta22] [alpha21; alpha22] [eta21; eta22] [mFuel21; mFuel22] [flapdeflection21; flapdeflection22] ThirdStageFlag];
dlmwrite(strcat('SecondStageStates',namelist{j}),['time (s) ' 'altitude (m) ' 'Mach ' 'longitude (rad) ' 'latitude (rad) ' 'velocity (m/s) ' 'trajectory angle (rad) ' 'heading angle (rad) ' 'angle of attack (rad) ' 'bank angle (rad) ' 'fuel mass (kg) ' 'Flap Deflection (deg)' 'Third Stage (flag)'],'');
dlmwrite(strcat('SecondStageStates',namelist{j}),SecondStageStates,'-append','delimiter',' ');
movefile(strcat('SecondStageStates',namelist{j}),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));




%% Check KKT and pontryagins minimum
% Check that the hamiltonian = 0 (for free end time)
% Necessary condition
input_test = output{j}.result.solution;
input_test.auxdata = auxdata;
phaseout_test = CombinedContinuous(input_test);

H1 = [];
H2 = [];
H3 = [];
H4 = [];

lambda1 = output{j}.result.solution.phase(1).costate;
for i = 1:length(lambda1)-1
    H1(i) = lambda1(i+1,:)*phaseout_test(1).dynamics(i,:).'; %H = lambda transpose * f(x,u,t) + L, note that there is no continuous cost L
end

lambda2 = output{j}.result.solution.phase(2).costate;
for i = 1:length(lambda2)-1
    H2(i) = lambda2(i+1,:)*phaseout_test(2).dynamics(i,:).'; %H = lambda transpose * f(x,u,t) + L, note that there is no continuous cost L
end

lambda3 = output{j}.result.solution.phase(3).costate;
for i = 1:length(lambda3)-1
    H3(i) = lambda3(i+1,:)*phaseout_test(3).dynamics(i,:).'; %H = lambda transpose * f(x,u,t) + L, note that there is no continuous cost L
end

lambda4 = output{j}.result.solution.phase(4).costate;
for i = 1:length(lambda4)-1
    H4(i) = lambda4(i+1,:)*phaseout_test(4).dynamics(i,:).'; %H = lambda transpose * f(x,u,t) + L, note that there is no continuous cost L
end

figure(2410)
% hold on
% plot(time1(1:end-1),H1)
% plot(time21(1:end-1),H2)
% plot(time22(1:end-1),H3)
% plot(time3(1:end-1),H4)
% ylabel('Hamiltonian')
% xlabel('Time (s)')
% legend('First Stage Stage Ascent','Second Stage Ascent','Second Stage Return','Third Stage Powered Ascent')




subplot(2,2,1)
hold on
title('First Stage')
plot(time1(1:end-1)-time1(1),H1);
ylabel('Hamiltonian')
xlabel('Time (s)')
ylim([-1 1]);
xlim([0 time1(end-1)-time1(1)]);

subplot(2,2,2)
hold on
title('Second Stage Ascent')
plot(time21(1:end-1)-time21(1),H2);
ylabel('Hamiltonian')
xlabel('Time (s)')
ylim([-1 1]);
xlim([0 time21(end-1)-time21(1)]);
subplot(2,2,3)
hold on
title('Second Stage Return')
plot(time22(1:end-1)-time22(1),H3);
ylabel('Hamiltonian')
xlabel('Time (s)')
ylim([-1 1]);
xlim([0 time22(end-1)-time22(1)]);
subplot(2,2,4)
hold on
title('Third Stage')
plot(time3(1:end-1)-time3(1),H4);
ylabel('Hamiltonian')
xlabel('Time (s)')
ylim([-1 1]);
xlim([0 time3(end-1)-time3(1)]);

%% Check State Feasibility
% Check calculated derivatives with the numerical derivative of each
% porimal, scaled by that primal

% for i = 1:length(output{j}.result.solution.phase(2).state(1,:))
% plot(time21(1:end-1),([diff(output{j}.result.solution.phase(2).state(:,i))./diff(output{j}.result.solution.phase(2).time)] - phaseout_test(2).dynamics(1:end-1,i))./max(abs(phaseout_test(2).dynamics(1:end-1,i))),'--');
% end
% for i = 1:length(output{j}.result.solution.phase(3).state(1,:))
%     if i<= 7 % Plot different line styles when no. of colours exceeded
%     plot(time22(1:end-1),([diff(output{j}.result.solution.phase(3).state(:,i))./diff(output{j}.result.solution.phase(3).time)] - phaseout_test(3).dynamics(1:end-1,i))./max(abs(phaseout_test(3).dynamics(1:end-1,i))));
%     else
%     plot(time22(1:end-1),([diff(output{j}.result.solution.phase(3).state(:,i))./diff(output{j}.result.solution.phase(3).time)] - phaseout_test(3).dynamics(1:end-1,i))./max(abs(phaseout_test(3).dynamics(1:end-1,i))),':');
%     end
% end

% 
% for i = 1:length(output{j}.result.solution.phase(2).state(1,:))
% plot(time21(1:end),(cumtrapz(output{j}.result.solution.phase(2).time,phaseout_test(2).dynamics(1:end,i)) + output{j}.result.solution.phase(2).state(1,i)- output{j}.result.solution.phase(2).state(:,i))./(max(output{j}.result.solution.phase(2).state(:,i))),'--');
% end
% for i = 1:length(output{j}.result.solution.phase(3).state(1,:))
%     if i<= 7 % Plot different line styles when no. of colours exceeded
%     plot(time22(1:end),(cumtrapz(output{j}.result.solution.phase(3).time,phaseout_test(3).dynamics(1:end,i)) + output{j}.result.solution.phase(3).state(1,i)- output{j}.result.solution.phase(3).state(:,i))./(max(output{j}.result.solution.phase(3).state(:,i))),'-');
%     else
%     plot(time22(1:end),(cumtrapz(output{j}.result.solution.phase(3).time,phaseout_test(3).dynamics(1:end,i)) + output{j}.result.solution.phase(3).state(1,i)- output{j}.result.solution.phase(3).state(:,i))./(max(output{j}.result.solution.phase(3).state(:,i))),':');
%     end
% end

figure(2420)
fig = gcf;
set(fig,'Position',[200 200 700 600]);
hold on

%suptitle('Dynamics Feasibility Check')

for i = [1:6 9]
    if i == 1
        i1 = 1;
    elseif i ==2
        i1 = 9;
   elseif i ==3
        i1 = 8;
   elseif i ==4
        i1 = 2;
   elseif i ==5
        i1 = 4;
   elseif i ==6
        i1 = 6;
   elseif i ==9
        i1 = 3;
    end
   
        
    Stage1_int = cumtrapz(output{j}.result.solution.phase(1).time,phaseout_test(1).dynamics(1:end,i1));
    Stage1_error = (Stage1_int+ output{j}.result.solution.phase(1).state(1,i1)- output{j}.result.solution.phase(1).state(:,i1))./(max(output{j}.result.solution.phase(1).state(:,i1))-min(output{j}.result.solution.phase(1).state(:,i1)));
%     Stage1_int = cumtrapz([time1(1):0.1:time1(end)],pchip(output{j}.result.solution.phase(1).time,phaseout_test(1).dynamics(1:end,i),[time1(1):0.1:time1(end)]))';
%     Stage1_error = (interp1([time1(1):0.1:time1(end)],Stage1_int,time1)'+ output{j}.result.solution.phase(1).state(1,i)- output{j}.result.solution.phase(1).state(:,i))./(max(output{j}.result.solution.phase(1).state(:,i))-min(output{j}.result.solution.phase(1).state(:,i)));

%     Stage2_int = cumtrapz(output{j}.result.solution.phase(2).time,phaseout_test(2).dynamics(1:end,i));
%    Stage2_error = (Stage2_int+ output{j}.result.solution.phase(2).state(1,i)- output{j}.result.solution.phase(2).state(:,i))./(max(output{j}.result.solution.phase(2).state(:,i))-min(output{j}.result.solution.phase(2).state(:,i)));
Stage2_int = cumtrapz([time21(1):0.1:time21(end)],pchip(output{j}.result.solution.phase(2).time,phaseout_test(2).dynamics(1:end,i),[time21(1):0.1:time21(end)]))';
    Stage2_error = (interp1([time21(1):0.1:time21(end)],Stage2_int,time21)+ output{j}.result.solution.phase(2).state(1,i)- output{j}.result.solution.phase(2).state(:,i))./(max(output{j}.result.solution.phase(2).state(:,i))-min(output{j}.result.solution.phase(2).state(:,i)));
%     Return_int = cumtrapz(output{j}.result.solution.phase(3).time,phaseout_test(3).dynamics(1:end,i));
%     Return_error = (Return_int+ output{j}.result.solution.phase(3).state(1,i)- output{j}.result.solution.phase(3).state(:,i))./(max(output{j}.result.solution.phase(3).state(:,i))-min(output{j}.result.solution.phase(3).state(:,i)));

if returnMode == 1
Return_int = cumtrapz([time22(1):0.1:time22(end)],pchip(output{j}.result.solution.phase(3).time,phaseout_test(3).dynamics(1:end,i),[time22(1):0.1:time22(end)]))';
    Return_error = (interp1([time22(1):0.1:time22(end)],Return_int,time22)+ output{j}.result.solution.phase(3).state(1,i)- output{j}.result.solution.phase(3).state(:,i))./(max(output{j}.result.solution.phase(3).state(:,i))-min(output{j}.result.solution.phase(3).state(:,i)));
else
 Return_error = 0;  
end
    
    
    if i == 1
        i3 = 1;
    elseif i ==2
        i3 = 0;
   elseif i ==3
        i3 = 6;
   elseif i ==4
        i3 = 2;
   elseif i ==5
        i3 = 3;
   elseif i ==6
        i3 = 7;
   elseif i ==9
        i3 = 4;
    end
    
    if i3 == 0
        Stage3_error = zeros(length(time3),1);
    else
    Stage3_int = cumtrapz(output{j}.result.solution.phase(4).time,phaseout_test(4).dynamics(1:end,i3));
    Stage3_error = (Stage3_int+ output{j}.result.solution.phase(4).state(1,i3)- output{j}.result.solution.phase(4).state(:,i3))./(max(output{j}.result.solution.phase(4).state(:,i3))-min(output{j}.result.solution.phase(4).state(:,i3)));
    
%     Stage3_int = cumtrapz([time3(1):0.1:time3(end)],pchip(output{j}.result.solution.phase(4).time,phaseout_test(4).dynamics(1:end,i),[time3(1):0.1:time3(end)]))';
%     Stage3_error = (interp1([time3(1):0.1:time3(end)],Stage3_int,time3)+ output{j}.result.solution.phase(4).state(1,i)- output{j}.result.solution.phase(4).state(:,i))./(max(output{j}.result.solution.phase(4).state(:,i))-min(output{j}.result.solution.phase(4).state(:,i)));

    end

    
% plot([time1'; time21 ;time22; time3],[Stage1_error;Stage2_error; Return_error;Stage3_error]);


subplot(2,2,1)
hold on
title('First Stage')
plot(time1',Stage1_error*100);
ylabel('Integrated Error (%)')
xlabel('Time (s)')
ylim([-1 1]);
xlim([0 time1(end)-time1(1)]);

subplot(2,2,2)
hold on
title('Second Stage Ascent')
plot(time21-time21(1),Stage2_error*100);
ylabel('Integrated Error (%)')
xlabel('Time (s)')
ylim([-1 1]);
xlim([0 time21(end)-time21(1)]);
subplot(2,2,3)
hold on
title('Second Stage Return')
plot(time22-time22(1),Return_error*100);
ylabel('Integrated Error (%)')
xlabel('Time (s)')
ylim([-1 1]);
xlim([0 time22(end)-time22(1)]);
subplot(2,2,4)
hold on
title('Third Stage')
plot(time3-time3(1),Stage3_error*100);
ylabel('Integrated Error (%)')
xlabel('Time (s)')
ylim([-1 1]);
xlim([0 time3(end)-time3(1)]);
end




% ylim([-1,1])
% legend('Alt Ascent','lon Ascent','lat Ascent','v Ascent','gamma Ascent','zeta Ascent','aoa Ascent','bank Ascent','mFuel Ascent', 'Alt Descent','lon Descent','lat Descent','v Descent','gamma Descent','zeta Descent','aoa Descent','bank Descent','mFuel Descent','throttle Descent')
legend('Alt ','lon ','lat','v ','gamma ','zeta ','Mass');


%% Plot Heating Rate

figure(701)
plot(time21,heating_rate21/1000);
hold on
plot(time22,heating_rate22/1000);
xlabel('time');
ylabel('Heat Flux (kW/m^2)');
title('Stagnation Point Heating')


%% plot engine interpolation visualiser
T0 = spline( auxdata.interp.Atmosphere(:,1),  auxdata.interp.Atmosphere(:,2), alt21); 
T_in1 = auxdata.interp.tempgridded(M21,rad2deg(alpha21)).*T0;
M_in1 = auxdata.interp.M1gridded(M21, rad2deg(alpha21));

plotM = [min(M_englist):0.01:9];
plotT = [min(T_englist):1:600];
[gridM,gridT] =  ndgrid(plotM,plotT);
interpeq = auxdata.interp.eqGridded(gridM,gridT);
interpIsp = auxdata.interp.IspGridded(gridM,gridT);

T02 = spline( auxdata.interp.Atmosphere(:,1),  auxdata.interp.Atmosphere(:,2), alt22); 
T2_in1 = auxdata.interp.tempgridded(M22,rad2deg(alpha22)).*T02;
M2_in1 = auxdata.interp.M1gridded(M22, rad2deg(alpha22));

IspPlot = 1
if IspPlot == 1

figure(2100)
hold on
contourf(gridM,gridT,interpeq,100,'LineColor','none');
scatter(engine_data(:,1),engine_data(:,2),30,engine_data(:,4),'k');
xlabel('M1')
ylabel('T1')
c=colorbar
c.Label.String = 'Equivalence Ratio';
caxis([.4 1])
plot(M_in1,T_in1,'r');

error_Isp = auxdata.interp.IspGridded(engine_data(:,1),engine_data(:,2))-engine_data(:,3);
% error_Isp =  polyvaln(auxdata.interp.Isppoly,[engine_data(:,1),engine_data(:,2)])-engine_data(:,3)


figure(2110)
hold on
contourf(gridM,gridT,interpIsp,100,'LineColor','none');
scatter(engine_data(:,1),engine_data(:,2),30,engine_data(:,3),'k')
xlabel('M1')
ylabel('T1')
c=colorbar
caxis([0 3000])
c.Label.String = 'ISP';
plot(M_in1,T_in1,'r');

figure(2210)
hold on
contourf(gridM,gridT,interpIsp,100,'LineColor','none');
scatter(engine_data(:,1),engine_data(:,2),30,engine_data(:,3),'k')
xlabel('M1')
ylabel('T1')
c=colorbar
caxis([0 3000])
c.Label.String = 'ISP';
plot(M2_in1,T2_in1,'r');
plot(M2_in1(throttle22>0.5 & time22>100+525 & time22<300+525),T2_in1(throttle22>0.5 & time22>100+525 & time22<300+525),'c');
plot(M2_in1(throttle22>0.1 & time22>300+525 & time22<450+525),T2_in1(throttle22>0.1 & time22>300+525 & time22<450+525),'c');
plot(M2_in1(throttle22>0.1 & time22>450+525),T2_in1(throttle22>0.1 & time22>450+525),'c');
figure(2120)
contourf(MList_EngineOn,AOAList_EngineOn,auxdata.interp.M1gridded(MList_EngineOn,AOAList_EngineOn),100,'LineWidth',0.01)
xlabel('M')
ylabel('Angle of Attack (deg)')
c=colorbar
c.Label.String = 'M1';

figure(2130)
contourf(MList_EngineOn,AOAList_EngineOn,auxdata.interp.tempgridded(MList_EngineOn,AOAList_EngineOn),100,'LineWidth',0.01)
xlabel('M')
ylabel('Angle of Attack (deg)')
c=colorbar
c.Label.String = 'T1/T0';

figure(2140)
contourf(MList_EngineOn,AOAList_EngineOn,auxdata.interp.presgridded(MList_EngineOn,AOAList_EngineOn),100,'LineWidth',0.01)
xlabel('M')
ylabel('Angle of Attack (deg)')
c=colorbar
c.Label.String = 'P1/P0';

end

%%
[gridM2,gridAoA2] =  ndgrid(plotM,plotT);


%% ThirdStage


forward0 = [alt3(1),v3(1),gamma3(1),m3(1),lat3(1),zeta3(1)];

% [f_t, f_y] = ode45(@(f_t,f_y) ForwardSim(f_y,AlphaInterp(t,Alpha,f_t),communicator,communicator_trim,SPARTAN_SCALE,Atmosphere,mode,scattered),t,forward0);
[f_t, f_y] = ode45(@(f_t,f_y) VehicleModel3_forward(f_t, f_y,auxdata,ControlInterp(time3,aoa3,f_t),ControlInterp(time3,aoadot3,f_t)),time3(1:end),forward0);

% [rdot3,xidot3,phidot3,gammadot3,vdot3,zetadot3, mdot3, Vec_angle3, AoA_max3, T3, L3, D3, q3] = ThirdStageDynamics(alt3,gamma3,v3,m3,aoa3,time3,auxdata,aoadot3,lat3,zeta3);

% lon3 = [];
% lon3(1) = lon21(end);
% for i = 2:length(time3)
%     lon3(i) = lon3(i-1) + xidot3(i-1)*(time3(i)-time3(i-1));
% end


% [AltF_actual, v3F, altexo, v3exo, timeexo, mpayload, Alpha3, mexo,qexo,gammaexo,Dexo,zetaexo,latexo,incexo,Texo,CLexo,Lexo,incdiffexo,lonexo] = ThirdStageSim(alt3(end),gamma3(end),v3(end), lat3(end),lon3(end), zeta3(end), m3(end), auxdata);


forward_error3 = [(alt3-f_y(:,1))/(max(alt3)-min(alt3)) zeros(length(alt3),0)...
    (lat3-f_y(:,5))/(max(lat3)-min(lat3)) (v3-f_y(:,2))/(max(v3)-min(v3))...
    (gamma3-f_y(:,3))/(max(gamma3)-min(gamma3)) (zeta3-f_y(:,6))/(max(zeta3)-min(zeta3)) (m3-f_y(:,4))/(max(m3)-min(m3))];



figure(312)
hold on
plot(f_t(1:end),f_y(:,1));
plot(time3,alt3);

figure(313)
hold on
plot(f_t(1:end),f_y(:,2));
plot(time3,v3);



% figure(311)
%     fig = gcf;
% set(fig,'Position',[200 0 850 650])
% %%suptitle('Third Stage Trajectory');
%     hold on
%     
%     subplot(4,2,1)
%     hold on
%     title('Altitude (km','FontSize',9);
%     plot([time3; timeexo.'+time3(end)], [alt3; altexo.']/1000,'Color','k')
%     xlim([time3(1) timeexo(end)+time3(end)])
% set(gca,'xticklabels',[])
%     subplot(4,2,2)
%     hold on
%     title('Dynamic Pressure (kPa','FontSize',9);
%     plot([time3; timeexo.'+time3(end)],[q3;qexo.';qexo(end)]/1000,'Color','k')
%     xlim([time3(1) timeexo(end)+time3(end)])
% set(gca,'xticklabels',[])
%     subplot(4,2,3)
%     hold on
%     title('Angle of Attack (deg)','FontSize',9);
%     plot([time3; timeexo.'+time3(end)],[rad2deg(aoa3);0*ones(length(timeexo),1)],'Color','k')
%     xlim([time3(1) timeexo(end)+time3(end)])
% set(gca,'xticklabels',[])
%     subplot(4,2,4)
%     hold on
%     title('Velocity (m/s)','FontSize',9);
%     plot([time3; timeexo.'+time3(end)],[v3;v3exo.'],'Color','k')
%     xlim([time3(1) timeexo(end)+time3(end)])
% set(gca,'xticklabels',[])
%     subplot(4,2,5)
%     hold on
%     title('Mass (kg)','FontSize',9);
%     plot([time3; timeexo.'+time3(end)],[ m3;mexo.';mexo(end)],'Color','k')
%     xlim([time3(1) timeexo(end)+time3(end)])
% set(gca,'xticklabels',[])
%     subplot(4,2,6)
%     hold on
%     title('Thrust Vector Angle (deg)','FontSize',9);
%     plot([time3; timeexo.'+time3(end)],[rad2deg(Vec_angle3);0*ones(length(timeexo),1)],'Color','k')
%     xlabel('Time (s)','FontSize',9);
%     xlim([time3(1) timeexo(end)+time3(end)])
% set(gca,'xticklabels',[])
%     subplot(4,2,7)
%     hold on
%     title('Thrust','FontSize',9);
%     plot([time3], [T3],'Color','k')
% 
%     xlabel('Time (s)','FontSize',9);
%     xlim([time3(1) timeexo(end)+time3(end)])
%     subplot(4,2,8)
%     hold on
%     title('Drag','FontSize',9);
%     plot([time3], [D3],'Color','k')
% 
%     xlabel('Time (s)','FontSize',9);
%     xlim([time3(1) timeexo(end)+time3(end)])
%     
    
    
 figure(311)
 fig = gcf;
set(fig,'Position',[200 0 850 720])
hold on
% %suptitle('Second Stage Ascent')
subplot(3,1,1)
% set(gca,'xticklabels',[])
hold on
xlim([0 timeexo(end).'+time3(end)-time3(1)]);
 plot([time3-time3(1); timeexo.'+time3(end)-time3(1)],[alt3; altexo.']/1000,'-','color','k', 'linewidth', 1.);
% ylim([-30,40]);
ylabel('altitude(km)');
addaxis([time3-time3(1) ; timeexo.'+time3(end)-time3(1)],[v3;v3exo.']/1000,'--','color','k', 'linewidth', 1.);
addaxislabel(2,'Velocity (km/s)');

% addaxis(time,fpa,':','color','k', 'linewidth', 1.);
% addaxislabel(3,'Trajectory Angle (deg)');


addaxis([time3-time3(1); timeexo.'+time3(end)-time3(1)],[q3;qexo.';qexo(end)]/1000,':','color','k', 'linewidth', 1.2);
addaxislabel(3,'Dynamic Pressure (kPa)');

legend(  'Altitude', 'Velocity', 'Dynamic Pressure', 'location', 'best');


subplot(3,1,2)
% set(gca,'xticklabels',[])
hold on
xlim([0 timeexo(end).'+time3(end)-time3(1)]);
 plot([time3-time3(1); timeexo.'+time3(end)-time3(1)],[rad2deg(gamma3);rad2deg(gammaexo)'],'-','color','k', 'linewidth', 1.);
% ylim([-30,40]);
ylabel('Trajectory Angle (deg)');

addaxis([time3-time3(1); timeexo.'+time3(end)-time3(1)],[rad2deg(zeta3);rad2deg(zetaexo)'],'--','color','k', 'linewidth', 1.);
addaxislabel(2,'Heading Angle (deg)');

addaxis([time3-time3(1); timeexo.'+time3(end)-time3(1)],[ m3;mexo.';mexo(end)],':','color','k', 'linewidth', 1.2);
addaxislabel(3,'Mass (kg)');

legend(  'Trajectory Angle', 'Heading Angle', 'Mass', 'location', 'best');


subplot(3,1,3)
% set(gca,'xticklabels',[])
hold on
xlim([0 timeexo(end).'+time3(end)-time3(1)]);
 plot([time3-time3(1); timeexo.'+time3(end)-time3(1)],[rad2deg(aoa3);0*ones(length(timeexo),1)],'-','color','k', 'linewidth', 1.);
% ylim([-30,40]);
ylabel('Angle of Attack (deg)');

addaxis([time3-time3(1); timeexo.'+time3(end)-time3(1)],[T3/1000;0*ones(length(timeexo),1)],'--','color','k', 'linewidth', 1.);
addaxislabel(2,'Thrust (kN)');

addaxis([time3-time3(1); timeexo.'+time3(end)-time3(1)],[rad2deg(Vec_angle3);0*ones(length(timeexo),1)],':','color','k', 'linewidth', 1.2);
addaxislabel(3,'Thrust Vector Angle (deg)');



legend(  'Angle of Attack', 'Thrust (kN)', 'Thrust Vector', 'location', 'best');

xlabel('Time (s)');


    
    % Write data to file
    dlmwrite(strcat('ThirdStageStates',namelist{j}),['time (s) ' 'altitude (m) ' 'velocity (m/s) ' 'mass (kg) ' 'dynamic pressure (Pa)' 'trajectory angle (rad) ' 'Lift (N)' 'Drag (N)' 'heading angle (rad) ' 'latitude (rad) ' 'angle of attack (rad) '],'');
    dlmwrite(strcat('ThirdStageStates',namelist{j}),[[time3; time3(end)+timeexo'], [alt3; altexo'], [v3; v3exo'], [m3; mexo'; mexo(end)],[q3; qexo'; qexo(end)] ,[gamma3; gammaexo'],[L3; Lexo'; Lexo(end)],[D3; Dexo'; Dexo(end)] ,[zeta3; zetaexo'], [lat3; latexo'], [aoa3; zeros(length(timeexo),1)]],'-append','delimiter',' ')
movefile(strcat('ThirdStageStates',namelist{j}),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));

% 
%% First Stage =========================================================

% 
FirstStageSMF = (mRocket - mFuel)/(m1(1) - mSpartan);
% 

FirstStageStates = [time1' alt1' v1' m1' gamma1' alpha1' zeta1' lat1' lon1'];

dlmwrite(strcat('FirstStageStates',namelist{j}),['time (s) ' 'altitude (m) ' 'velocity (m/s) ' 'mass (kg)' 'trajectory angle (rad) ' 'angle of attack (rad) ' 'heading angle (rad) ' 'latitude (rad)'  'longitude (rad)'],'');
dlmwrite(strcat('FirstStageStates',namelist{j}),FirstStageStates,'-append','delimiter',' ');
movefile(strcat('FirstStageStates',namelist{j}),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));


% Iterative Prepitch Determination ========================================
%This back determines the mass and launch altitude necessary to get to
%100m, 30m/s at the PS method determined fuel mass

interp = auxdata.interp;
% Throttle = auxdata.Throttle;
Vehicle = auxdata.Vehicle;
Atmosphere = auxdata.Atmosphere;

% ntoe that launch altitude does vary, but it should only be slightly
% controls = fminunc(@(controls) prepitch(controls,m1(1),interp,Throttle,Vehicle,Atmosphere,auxdata),[-50,10]);
controls = fminunc(@(controls) prepitch(controls,m1(1),interp,Vehicle,Atmosphere,auxdata),[-50,10]);

h_launch = controls(1)
t_prepitch = controls(2)
Isp1 = Vehicle.Isp.SL;
T121 = Vehicle.T.SL.*1; % Throttle is set to 1 for prepitch
dm1 = -T121./Isp1./9.81;
m0_prepitch = m1(1) - dm1*t_prepitch;



%% Forward Integrator
 phase = 'postpitch';
tspan = time1; 
% postpitch0_f = [y(end,1) y(end,2) y(end,3) deg2rad(89.9) phi(1) zeta(1)]; % set mass
postpitch0_f = [h0 v0 m1(1) deg2rad(89.9) lat1(1) zeta1(1)];

% [t_postpitch_f, postpitch_f] = ode45(@(t_f,postpitch_f) FirstStageDynamicsForward(postpitch_f,ControlFunction(t_f,time1,zeta1),ControlFunction(t_f,time1,alpha1),phase,interp,Throttle,Vehicle,Atmosphere,auxdata), tspan, postpitch0_f);
% [t_postpitch_f, postpitch_f] = ode45(@(t_f,postpitch_f) FirstStageDynamicsForward(postpitch_f,ControlFunction(t_f,time1,zeta1),ControlFunction(t_f,time1,alpha1),phase,interp,Vehicle,Atmosphere,auxdata), tspan, postpitch0_f);
[t_postpitch_f, postpitch_f] = ode45(@(t_f,postpitch_f) FirstStageDynamicsForward(postpitch_f,ControlFunction(t_f,time1,zeta1),ControlFunction(t_f,time1,alpha1),phase,interp,ControlFunction(t_f,time1,Throttle1),Vehicle,Atmosphere,auxdata), tspan, postpitch0_f);


forward_error1 = [(alt1'-postpitch_f(:,1))/(max(alt1)-min(alt1)) zeros(length(alt1),1)...
    (lat1'-postpitch_f(:,5))/(max(lat1)-min(lat1)) (v1'-postpitch_f(:,2))/(max(v1)-min(v1))...
    (gamma1'-postpitch_f(:,4))/(max(gamma1)-min(gamma1)) (zeta1'-postpitch_f(:,6))/(max(zeta1)-min(zeta1)) (m1'-postpitch_f(:,3))/(max(m1)-min(m1))];


figure(103)
hold on
plot(postpitch_f(:,1));
plot(alt1);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Pre-Pitchover Simulation                         %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
h0_prepitch = h_launch;  %Rocket starts on the ground
v0_prepitch = 0;  %Rocket starts stationary
gamma0_prepitch = deg2rad(90);

phase = 'prepitch';
tspan2 = [0 t_prepitch]; % time to fly before pitchover (ie. straight up)

y0 = [h0_prepitch, v0_prepitch, m0_prepitch, gamma0_prepitch, 0, 0, 0, 0, 0, 1]

% this performs a forward simulation before pitchover. The end results of
% this are used as initial conditions for the optimiser. 
% [t_prepitch, y] = ode45(@(t,y) FirstStageDynamics(y,0,0,phase,interp,Throttle,Vehicle,Atmosphere,auxdata), tspan2, y0);  
[t_prepitch, y] = ode45(@(t,y) FirstStageDynamics(y,[0; 1],0,phase,interp,Vehicle,Atmosphere,auxdata), tspan2, y0);  


% figure(111);
% hold on
% %%suptitle('First Stage Trajectory');
%     fig = gcf;
% set(fig,'Position',[200 0 850 600])
% subplot(4,2,1)
% hold on
% title('Trajectory Angle (deg)');
% xlim([0 time1(end)+t_prepitch(end)]);
% plot([t_prepitch.' time1+t_prepitch(end)], [rad2deg(y(:,4).') rad2deg(gamma1)],'color','k');
% subplot(4,2,2)
% hold on
% title('Velocity (m/s)');
% xlim([0 time1(end)+t_prepitch(end)]);
% plot([t_prepitch.' time1+t_prepitch(end)], [y(:,2).' v1],'color','k');
% subplot(4,2,3)
% hold on
% title('Altitude (km)');
% xlim([0 time1(end)+t_prepitch(end)]);
% plot([t_prepitch.' time1+t_prepitch(end)], [y(:,1).'/1000 alt1/1000],'color','k');
% subplot(4,2,4)
% hold on
% title('Angle of Attack (deg)');
% xlim([0 time1(end)+t_prepitch(end)]);
% plot([t_prepitch.' time1+t_prepitch(end)], [zeros(1,length(t_prepitch)) rad2deg(alpha1)],'color','k');
% subplot(4,2,5)
% hold on
% title('Mass (kg)');
% xlim([0 time1(end)+t_prepitch(end)]);
% plot([t_prepitch.' time1+t_prepitch(end)], [y(:,3).' m1],'color','k');
% subplot(4,2,6)
% hold on
% title('Heading Angle (deg)');
% xlim([0 time1(end)+t_prepitch(end)]);
% plot([time1+t_prepitch(end)], [rad2deg(zeta1)],'color','k');
% subplot(4,2,7)
% hold on
% title('Latitude (deg)');
% xlim([0 time1(end)+t_prepitch(end)]);
% plot([t_prepitch.' time1+t_prepitch(end)], [rad2deg(lat1(1)+y(:,8).') rad2deg(lat1)],'color','k');
% 
% % plot([primal.nodes], [rad2deg(gamma)/100],'color','k','linestyle','-');
% % plot([primal.nodes], [v/1000],'color','k','linestyle','--');
% % plot([primal.nodes], [V/10000],'color','k','linestyle',':');
% % plot([primal.nodes], [rad2deg(alpha)/10],'color','k','linestyle','-.')
% xlabel('Time (s)')
% xlim([0,time1(end)+t_prepitch(end)]);



phase = 'postpitch';

interp = auxdata.interp;
% Throttle = auxdata.Throttle;
Vehicle = auxdata.Vehicle;
Atmosphere = auxdata.Atmosphere;

% [dz1,q1,xi1,vec_angle1,T1,D1,dm1] = FirstStageDynamics(output{j}.result.solution.phase(1).state',output{j}.result.solution.phase(1).control',output{j}.result.solution.phase(1).time',phase,interp,Throttle,Vehicle,Atmosphere,auxdata); % Pass primal variables into dynamics simulator
[dz1,q1,xi1,vec_angle1,T121,D1,dm1,M1,rho01,P01,T01,P_1_tip1,T_1_tip1,rho_1_tip1,M_1_tip1,v_1_tip1,heating_rate1,heating_rate_LE1] = FirstStageDynamics(output{j}.result.solution.phase(1).state',output{j}.result.solution.phase(1).control',output{j}.result.solution.phase(1).time',phase,interp,Vehicle,Atmosphere,auxdata); % Pass primal variables into dynamics simulator


 figure(111)
 fig = gcf;
set(fig,'Position',[200 0 850 720])
hold on
% %suptitle('Second Stage Ascent')
subplot(3,1,1)
set(gca,'xticklabels',[])
hold on
xlim([0 time1(end)-time1(1)]);
 plot(time1-time1(1),alt1/1000,'-','color','k', 'linewidth', 1.);
% ylim([-30,40]);
ylabel('altitude(km)');
addaxis(time1-time1(1),v1/1000,'--','color','k', 'linewidth', 1.);
addaxislabel(2,'Velocity (km/s)');

% addaxis(time,fpa,':','color','k', 'linewidth', 1.);
% addaxislabel(3,'Trajectory Angle (deg)');


addaxis(time1-time1(1),q1/1000,':','color','k', 'linewidth', 1.2);
addaxislabel(3,'Dynamic Pressure (kPa)');

legend(  'Altitude', 'Velocity', 'Dynamic Pressure', 'location', 'best');


subplot(3,1,2)
set(gca,'xticklabels',[])
hold on
xlim([0 time1(end)-time1(1)]);
 plot(time1-time1(1),rad2deg(gamma1),'-','color','k', 'linewidth', 1.);
% ylim([-30,40]);
ylabel('Trajectory Angle (deg)');

addaxis(time1-time1(1),T121/1000,'--','color','k', 'linewidth', 1.);
ylabel('Thrust (kN)');

% xlim([0 time1(end)-time1(1)]);
% addaxis(time1-time1(1),rad2deg(alpha1),'-','color','k', 'linewidth', 1.);
% % ylim([-30,40]);
% addaxislabel(3,'Angle of Attack (deg)');


addaxis(time1-time1(1),rad2deg(alpha1),':','color','k', 'linewidth', 1.2);
addaxislabel(3,'Angle of Attack (deg)');

% xlabel('Time (s)');

legend(  'Trajectory Angle', 'Thrust', 'Angle of Attack', 'location', 'best');



subplot(3,1,3)

hold on
xlim([0 time1(end)-time1(1)]);
 plot(time1-time1(1),Throttle1,'-','color','k', 'linewidth', 1.);

ylabel('Throttle');

addaxis(time1-time1(1),(T121-D1)./-dm1,'--','color','k', 'linewidth', 1.);
addaxislabel(2,'Net Isp (s)');



addaxis(time1-time1(1),rad2deg(vec_angle1),':','color','k', 'linewidth', 1.2);
addaxislabel(3,'Thrust Vector Angle (deg)');

xlabel('Time (s)');

legend(  'Throttle', 'Net Isp', 'Thrust Vector Angle', 'location', 'best');








saveas(figure(111),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('FirstStage',namelist{j},'.fig')]);
print(figure(111),strcat('FirstStage',namelist{j}),'-dpng');
movefile(strcat('FirstStage',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));


%%

FlowProps = [[time1'; time21; time22] [alt1'; alt21; alt22] [M1'; M21;M22] [T01';T021;T022] [P01';P021;P022] [rho01';rho021;rho022]  [M_1_tip1';M_1_tip21';M_1_tip22'] [T_1_tip1';T_1_tip21';T_1_tip22'] [P_1_tip1';P_1_tip21';P_1_tip22'] [rho_1_tip1';rho_1_tip21';rho_1_tip22'] [v_1_tip1';v_1_tip21';v_1_tip22'] [heating_rate1';heating_rate21;heating_rate22] [heating_rate_LE1';heating_rate_LE21;heating_rate_LE22]  [ones(length(alt1),1);2*ones(length(alt21),1);3*ones(length(alt22),1)]];
dlmwrite(strcat('FlowProps',namelist{j}),['time (s) ' 'altitude (m) ' 'Mach ' 'Incoming Air Temp (K)' 'Incoming Air Pressure (Pa)' 'Incoming Air density (kg/m^3)' 'Mach no. After Normal Shock' 'Temperature After Normal Shock (K)'  'Pressure After Normal Shock (Pa)' 'Density After Normal Shock (kg/m^3)' 'Velocity After Normal Shock (m/s)' 'heating rate (w/m^2)' 'heating rate LE (w/m^2)' ' Phase (first-second-return)'],'');   
dlmwrite(strcat('FlowProps',namelist{j}),FlowProps,'-append','delimiter',' ');
movefile(strcat('FlowProps',namelist{j}),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));


%% Plot Forward Sim Error
figure(550)
fig = gcf;
set(fig,'Position',[200 200 700 600]);
%suptitle('Normalised Forward Simulation Error')
subplot(2,2,1)
hold on
for i = 1:length(forward_error1(1,:))
    plot(time1-time1(1),forward_error1(:,i)*100)
end
ylabel('Integrated Error (%)')
xlabel('Time (s)')
xlim([0 time1(end)-time1(1)]);
title('First Stage')

subplot(2,2,2)
hold on
for i = 1:length(forward_error21(1,:))
    plot(time21-time21(1),forward_error21(:,i)*100)
end
ylabel('Integrated Error (%)')
xlabel('Time (s)')
xlim([0 time21(end)-time21(1)]);
title('Second Stage Ascent')

subplot(2,2,3)
hold on
for i = 1:length(forward_error22(1,:))
    plot(time22-time22(1),forward_error22(:,i)*100)
end
ylabel('Integrated Error (%)')
xlabel('Time (s)')
xlim([0 time22(end)-time22(1)]);
title('Second Stage Return')
legend('Alt ','lon ','lat','v ','gamma ','zeta ','Mass','Location','best');

subplot(2,2,4)
hold on
for i = 1:length(forward_error3(1,:))
    plot(time3-time3(1),forward_error3(:,i)*100)
end
ylabel('Integrated Error (%)')
xlabel('Time (s)')
xlim([0 time3(end)-time3(1)]);
title('Third Stage Powered Ascent')

% hSub = subplot(2,3,3); plot(1, nan, 1, nan, 1, nan, 1, nan, 1, nan, 1, nan, 'r'); set(hSub, 'Visible', 'off');

%% Calculate Change in Exergy 
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7116765
% https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20180004205.pdf

H_RP1 = 43.2*10^6; %MJ/kg, from http://webserver.dmt.upm.es/~isidoro/bk3/c15/Fuel%20properties.pdf
rat_RP1 = 1/3.56;
H_H = 119.96*10^6; %https://h2tools.org/hyarc/calculator-tools/lower-and-higher-heating-values-fuels

dExergy_1 = ( (auxdata.Stage2.mStruct + mFuel21(1) +auxdata.Stage3.mTot)*v1(end)^2/2 - (auxdata.Stage2.mStruct + mFuel21(1) +auxdata.Stage3.mTot)*v1(1)^2/2) + cumtrapz((alt1+6370000),(auxdata.Stage2.mStruct + mFuel21(1) +auxdata.Stage3.mTot)*9.81.*(6370000./(alt1+6370000)).^2);

Spec_Exergy_1 = dExergy_1(end)/(m1(1)-m1(end));
Exergy_tot1 = (H_RP1*rat_RP1)*(m1(1)-m1(end));
eff_Exergy_1 = 1 - (Exergy_tot1 - dExergy_1)/(Exergy_tot1);

% Work done accelerating first stage structural mass
W_1= (m1(end) - auxdata.Stage2.mStruct - mFuel21(1) - auxdata.Stage3.mTot)*v1(end)^2/2 - (m1(end) - auxdata.Stage2.mStruct - mFuel21(1) - auxdata.Stage3.mTot)*v1(1)^2/2 + cumtrapz((alt1+6370000),(m1(end) - auxdata.Stage2.mStruct - mFuel21(1) - auxdata.Stage3.mTot)*9.81.*(6370000./(alt1+6370000)).^2);
W_1_pc = W_1/(Exergy_tot1)*100;

W_D1 = cumtrapz(time1,v1.*D1); % work done to overcome drag
W_D1_pc  = W_D1(end)/(Exergy_tot1)*100;
% work to accelerate fuel mass
P_1 = T121.*v1;% propulsive power
% W_P1 = cumtrapz(time1,P_1); % Propulsive work
% P_loss1 = Exergy_tot1 - W_P1 ; % losses due to propulsive inefficiencies
% P_loss1_pc  = P_loss1(end)/(Exergy_tot1)*100; %this seems to not quite add to 100%,probably due to component losses due to aoa and thrust vectoring
mF1 = m1 - m1(end);
W_mF1 = cumtrapz(v1,mF1.*v1) + cumtrapz((alt1+6370000),(mF1)*9.81.*(6370000./(alt1+6370000)).^2); % work done accelerating and lifting fuel mass
W_mF1_pc  = W_mF1(end)/(Exergy_tot1)*100;

P_loss1_pc  = 100 - eff_Exergy_1*100 - W_mF1_pc - W_D1_pc - W_1_pc;

% work done accelerating payload
W_1Payload= (mpayload)*v1(end)^2/2 - (mpayload)*v1(1)^2/2 + cumtrapz((alt1+6370000),(mpayload)*9.81.*(6370000./(alt1+6370000)).^2);
W_1Payload_pc = W_1Payload/(Exergy_tot1)*100;


W_1nextstage = ( (auxdata.Stage2.mStruct + mFuel21(1) +auxdata.Stage3.mTot - mpayload)*v1(end)^2/2 - (auxdata.Stage2.mStruct + mFuel21(1) +auxdata.Stage3.mTot - mpayload)*v1(1)^2/2) + cumtrapz((alt1+6370000),(auxdata.Stage2.mStruct + mFuel21(1) +auxdata.Stage3.mTot - mpayload)*9.81.*(6370000./(alt1+6370000)).^2);
W_1nextstage_pc = W_1nextstage/(Exergy_tot1)*100;

dExergy_21 = ((auxdata.Stage3.mTot)*v21(end)^2/2 - (auxdata.Stage3.mTot)*v21(1)^2/2) + cumtrapz((alt21+6370000),(auxdata.Stage3.mTot)*9.81.*(6370000./(alt21+6370000)).^2);
Spec_Exergy_21 = dExergy_21(end)/(mFuel21(1)-mFuel21(end));
if auxdata.returnMode == 0
Exergy_tot21 = (H_H)*(mFuel21(1)-0); % all fuel is used in exergy calculation
else
  Exergy_tot21 = (H_H)*(mFuel21(1)-mFuel21(end)); 
end


eff_Exergy_21 = 1 - (Exergy_tot21 - dExergy_21)/(Exergy_tot21);

Exergy_tot2 = (H_H)*(mFuel21(1)-0); 
eff_Exergy_2 = 1 - (Exergy_tot2 - dExergy_21)/(Exergy_tot2);

dExergy_21_3 = ((auxdata.Stage3.mTot)*v21(end)^2/2 - (auxdata.Stage3.mTot)*v21(1)^2/2) + cumtrapz((alt21+6370000),(auxdata.Stage3.mTot)*9.81.*(6370000./(alt21+6370000)).^2);

% Work done on the SPARTAN, lost at stage separation
W_2= (auxdata.Stage2.mStruct)*v21(end)^2/2 - (auxdata.Stage2.mStruct)*v21(1)^2/2 + cumtrapz((alt21+6370000),(auxdata.Stage2.mStruct)*9.81.*(6370000./(alt21+6370000)).^2);

W_2_pc = W_2/(Exergy_tot21)*100;

W_D21 = cumtrapz(time21,v21.*(Fd21-auxdata.T_spline_Rear(M21,rad2deg(alpha21),alt21/1000)*auxdata.A.*q21)); % work done to overcome drag
W_D21_pc  = W_D21(end)/(Exergy_tot21)*100;
% work to accelerate fuel mass
% P_21 = T21.*v21;% propulsive power % this seems to not quite add to 100%,
% probably due to component losses due to aoa and thrust vectoring
% probably due to 
% W_P21 = cumtrapz(time21,P_21); % Propulsive work
% P_loss21 = Exergy_tot21 - W_P21 ; % losses due to propulsive inefficiencies

W_mF21 = cumtrapz(v21,mFuel21.*v21) + cumtrapz((alt21+6370000),(mFuel21)*9.81.*(6370000./(alt21+6370000)).^2); % work done accelerating and lifting fuel mass
W_mF21_pc  = W_mF21(end)/(Exergy_tot21)*100;

P_loss21_pc  = 100 - eff_Exergy_21*100 - W_mF21_pc - W_D21_pc - W_2_pc; % work lost due to inefficiencies in the porpulsion system (remainder of work)

% work done accelerating payload
W_21Payload= (mpayload)*v21(end)^2/2 - (mpayload)*v21(1)^2/2 + cumtrapz((alt21+6370000),(mpayload)*9.81.*(6370000./(alt21+6370000)).^2);
W_21Payload_pc = W_21Payload/(Exergy_tot21)*100;

W_21nextstage = ( (auxdata.Stage3.mTot - mpayload)*v21(end)^2/2 - (auxdata.Stage3.mTot - mpayload)*v21(1)^2/2) + cumtrapz((alt21+6370000),(auxdata.Stage3.mTot - mpayload)*9.81.*(6370000./(alt21+6370000)).^2);
W_21nextstage_pc = W_21nextstage/(Exergy_tot21)*100;


if auxdata.returnMode == 1
    
      Exergy_tot22 = (H_H)*(mFuel22(1)-0); 
    
% Work done on the SPARTAN, lost at stage separation
W_22= (auxdata.Stage2.mStruct)*v22(end)^2/2 - (auxdata.Stage2.mStruct)*v22(1)^2/2 + cumtrapz((alt22+6370000),(auxdata.Stage2.mStruct)*9.81.*(6370000./(alt22+6370000)).^2);

W_22_pc = W_22/(Exergy_tot22)*100;

W_D22 = cumtrapz(time22,v22.*(Fd22-auxdata.T_spline_Rear(M22,rad2deg(alpha22),alt22/1000)*auxdata.A.*q22)); % work done to overcome drag
W_D22_pc  = W_D22(end)/(Exergy_tot22)*100;

W_mF22 = cumtrapz(v22,mFuel22.*v22) + cumtrapz((alt22+6370000),(mFuel22)*9.81.*(6370000./(alt22+6370000)).^2); % work done accelerating and lifting fuel mass
W_mF22_pc  = W_mF22(end)/(Exergy_tot22)*100;

P_loss22_pc  = 100 - W_mF22_pc - W_D22_pc - W_22_pc; % work lost due to inefficiencies in the porpulsion system (remainder of work)


end

% Energy_3 = mexo(end)*v3exo(end)^2/2 + ((altexo(end)+6370000)*mexo(end)*9.81.*(altexo(end)+6370000)/6370000);
dExergy_3_atm = mpayload*v3exo(1)^2/2-mpayload*v3(1)^2/2 + cumtrapz(alt3+6370000,mpayload*9.81.*(6370000./([alt3]+6370000)).^2); % Useful energy
% Add the energy change over the circularisation and hohmann transfer, in
% exoatmospheric cooridnate system
dExergy_3_exo = mpayload*(v3exo_coordchange+dvtot3)^2/2 - mpayload*v3exo_coordchange^2/2 - ((altexo(end))*mpayload*9.81.*(6370000./(altexo(end)+6370000)).^2) + ((566890)*mpayload*9.81.*(6370000./(566890+6370000)).^2);

dExergy_3 = dExergy_3_atm(end) + dExergy_3_exo;
%only use heating value for RP-1 (fuel)
Exergy_tot3 = (H_RP1*rat_RP1)*(auxdata.Stage3.mTot-m3_4); % Total available energy from fuel
eff_Exergy_3 = 1 - (Exergy_tot3 - dExergy_3)/(Exergy_tot3);

% energy used to accelerate and lift third stage structural mass
W_3_atm = (m3_4-mpayload)*v3exo(end)^2/2-(m3_4-mpayload)*v3(1)^2/2 + cumtrapz([alt3+6370000; altexo'+6370000],(m3_4-mpayload)*9.81.*(6370000./([alt3+6370000; altexo'+6370000]).^2)); % work done accelerating the third stage structural mass

W_3_exo = (m3_4-mpayload)*(v3exo_coordchange+dvtot3)^2/2 - (m3_4-mpayload)*v3exo_coordchange^2/2 - ((altexo(end))*(m3_4-mpayload)*9.81.*(6370000./(altexo(end)+6370000)).^2) + ((566890)*(m3_4-mpayload)*9.81.*(6370000./(566890+6370000)).^2);

W_3_pc = (W_3_atm+W_3_exo)/(Exergy_tot3)*100;

W_D3 = cumtrapz(time3,v3.*D3); % work done to overcome drag
W_D3_pc  = W_D3(end)/(Exergy_tot3)*100;
% work to accelerate fuel mass
% P_3 = T3.*v3;% propulsive power
% W_P3 = cumtrapz(time3,P_3); % Propulsive work
% P_loss3 = Exergy_tot3 - W_P3(end) ; % losses due to propulsive inefficiencies
% P_loss3_pc  = P_loss3(end)/(Exergy_tot3)*100*100;


W_mF3 = cumtrapz(v3,m3-m3_4-hs.mHS.*v3) + cumtrapz((alt3+6370000),(m3-m3_4-hs.mHS)*9.81.*(6370000./(alt3+6370000)).^2); % work done accelerating and lifting fuel mass (before heat shield)
W_mF3_pc  = W_mF3(end)/(Exergy_tot3)*100;

HT_loss_pc =  ((H_RP1*rat_RP1)*(m3(end)-m3_4) - dExergy_3_exo -W_3_exo/(Exergy_tot3)*100)/Exergy_tot3*100;% energy lost during hohmann transfer

W_HS3_pc = (hs.mHS.*hs.v^2 - hs.mHS.*v3(1)^2 - (alt3(1)+6370000)*(m3_4-mpayload)*9.81.*(6370000./(altexo(end)+6370000)).^2 + (hs.alt+6370000)*(m3_4-mpayload)*9.81.*(6370000./(hs.alt+6370000)).^2)/(Exergy_tot3)*100; %work needed to accelerate and lift heat shield
P_loss3_pc  = 100 - dExergy_3/(Exergy_tot3)*100 - W_mF3_pc - W_D3_pc - HT_loss_pc - W_HS3_pc - W_3_pc ;  % in atmosphere propulsive losses

%Calculate energy wasted at each staging manoeuvre
% stage_dEx1 = (m1(end)-auxdata.Stage2.mStruct-mFuel21(1)-auxdata.Stage3.mTot)*v1(end)^2/2 + alt1(end)*(m1(end)-auxdata.Stage2.mStruct-mFuel21(1)-auxdata.Stage3.mTot)*9.81.*(6370000./(alt1(end)+6370000)).^2;
% stage_dEx21 = (mFuel21(end)+auxdata.Stage2.mStruct)*v21(end)^2/2 + alt21(end)*(mFuel21(end)+auxdata.Stage2.mStruct)*9.81.*(6370000./(alt21(end)+6370000)).^2;
% stage_dEx3 = ((auxdata.Stage3.mTot - hs.mHS)*0.09)*(v3exo_coordchange+dvtot3)^2/2 + ((566890)*((auxdata.Stage3.mTot - hs.mHS)*0.09)*9.81.*(6370000./(566890+6370000)).^2) ; %heat shield is included in main dExergy
% eff_Exergy_total = 1-(Exergy_tot1 - dExergy_1(end) + stage_dEx1 + Exergy_tot21 - dExergy_21(end) + stage_dEx21 + Exergy_tot3 - dExergy_3 + stage_dEx3)/(Exergy_tot1+Exergy_tot21+Exergy_tot3);
% 
% eff_Exergy_2_3 = 1-(Exergy_tot21 - dExergy_21(end) + stage_dEx21 + Exergy_tot3(end) - dExergy_3 + stage_dEx3)/(Exergy_tot21+Exergy_tot3);

eff_Exergy_total = ((m3_4 - (auxdata.Stage3.mTot - hs.mHS)*0.09)*(v3exo_coordchange+dvtot3)^2/2 + (m3_4 - (auxdata.Stage3.mTot - hs.mHS)*0.09)*9.81*566890)/(Exergy_tot1+Exergy_tot21+Exergy_tot3)

% eff_Exergy_2_3 = 1-(Exergy_tot21 - dExergy_21(end) + stage_dEx21 + Exergy_tot3(end) - dExergy_3 + stage_dEx3)/(Exergy_tot21+Exergy_tot3);
eff_Exergy_2_3 =0;



%% calculate ground distance during return
dist_cov22= [];
for i = 1: length(lat22)-1
latlon1=[rad2deg(lat22(i)) rad2deg(lon22(i))];
latlon2=[rad2deg(lat22(i+1)) rad2deg(lon22(i+1))]; 
 dist_cov22(i) = lldistkm(latlon1,latlon2);

end
total_dist_cov22 = sum(dist_cov22);

dist_cov21= [];
for i = 1: length(lat21)-1
latlon1=[rad2deg(lat21(i)) rad2deg(lon21(i))];
latlon2=[rad2deg(lat21(i+1)) rad2deg(lon21(i+1))]; 
 dist_cov21(i) = lldistkm(latlon1,latlon2);

end
total_dist_cov21 = sum(dist_cov21);

dist_cov1 = [];
for i = 1: length(lat1)-1
latlon1=[rad2deg(lat1(i)) rad2deg(lon1(i))];
latlon2=[rad2deg(lat1(i+1)) rad2deg(lon1(i+1))]; 
 dist_cov1(i) = lldistkm(latlon1,latlon2);

end
total_dist_cov1 = sum(dist_cov1);


%% Create Easy Latex Inputs
if auxdata.returnMode == 1
    returnparam = {''};
elseif auxdata.mode == 1 || auxdata.mode == 90 
    returnparam = {'NoReturn'};
else
    returnparam = {''};
end

dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\PayloadToOrbit', namelist{j},returnparam ,'}{ ', num2str(round(ThirdStagePayloadMass,1),'%.1f') , '}'), '-append' , 'delimiter','','newline', 'pc')

if length(output)>1
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\PayloadVar', namelist{j},returnparam ,'}{ ', num2str(round((ThirdStagePayloadMass + output{3}.result.objective)/-output{3}.result.objective*100,2),'%.2f') , '}'), '-append' , 'delimiter','','newline', 'pc')
end

dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\totalExergyEff', namelist{j},returnparam ,'}{ ', num2str(eff_Exergy_total*100,'%.3f') , '}'), '-append' , 'delimiter','','newline', 'pc')
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\ExergyEffsecondandthird', namelist{j},returnparam ,'}{ ', num2str(eff_Exergy_2_3*100,'%.3f') , '}'), '-append' , 'delimiter','','newline', 'pc')


dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\firstsecondSeparationAlt', namelist{j},returnparam ,'}{ ', num2str(round(alt21(1)/1000,2),'%.2f') , '}'), '-append','delimiter','','newline', 'pc')
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\firstsecondSeparationv', namelist{j},returnparam ,'}{ ', num2str(round(v21(1),0)) , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\firstsecondSeparationgamma', namelist{j},returnparam ,'}{ ', num2str(rad2deg(gamma21(1)),'%.1f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\firstdExergy', namelist{j},returnparam ,'}{ ', num2str(dExergy_1(end)/10^9,'%.2f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\firstSpecExergy', namelist{j},returnparam ,'}{ ', num2str(Spec_Exergy_1/10^6,'%.2f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\firstExergyEff', namelist{j},returnparam ,'}{ ', num2str(eff_Exergy_1(end)*100,'%.3f') , '}'), '-append','delimiter','','newline', 'pc');

dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\firstWnextStage', namelist{j},returnparam ,'}{ ', num2str(W_1nextstage_pc(end) ,'%.3f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\firstWpayload', namelist{j},returnparam ,'}{ ', num2str(W_1Payload_pc(end) ,'%.3f') , '}'), '-append','delimiter','','newline', 'pc');


dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\firstExergyDestroyed', namelist{j},returnparam ,'}{ ', num2str((Exergy_tot1-dExergy_1(end))/10^9,'%.3f') , '}'), '-append','delimiter','','newline', 'pc');


dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\firstEnergy', namelist{j},returnparam ,'}{ ', num2str((Exergy_tot1)/10^9,'%.2f') , '}'), '-append','delimiter','','newline', 'pc');

dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\WDone', namelist{j},returnparam ,'}{ ', num2str(W_D1_pc(end) ,'%.2f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\Plossone', namelist{j},returnparam ,'}{ ', num2str(P_loss1_pc(end) ,'%.2f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\WmFone', namelist{j},returnparam ,'}{ ', num2str(W_mF1_pc(end) ,'%.2f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\Wone', namelist{j},returnparam ,'}{ ', num2str(W_1_pc(end) ,'%.2f') , '}'), '-append','delimiter','','newline', 'pc');

dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\PlossoneCombined', namelist{j},returnparam ,'}{ ', num2str(P_loss1_pc(end)+W_mF1_pc(end) ,'%.2f') , '}'), '-append','delimiter','','newline', 'pc');


dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\FirstStageDist', namelist{j},returnparam ,'}{ ', num2str(round(total_dist_cov1,1),'%.1f') , '}'), '-append','delimiter','','newline', 'pc');




dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\FirstStageSMF', namelist{j},returnparam ,'}{ ', num2str(round(FirstStageSMF,3),'%.3f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\secondthirdSeparationAlt', namelist{j},returnparam ,'}{ ', num2str(round(alt21(end)/1000,2),'%.2f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\secondthirdSeparationv', namelist{j},returnparam ,'}{ ', num2str(round(v21(end),0)) , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\secondthirdSeparationq', namelist{j},returnparam ,'}{ ', num2str(round(q21(end)/1000,2),'%.1f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\secondthirdSeparationgamma', namelist{j},returnparam ,'}{ ', num2str(rad2deg(gamma21(end)),'%.1f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\secondthirdSeparationLD', namelist{j},returnparam ,'}{ ', num2str(round(LD21(end),2),'%.1f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\secondFlightTime', namelist{j},returnparam ,'}{ ', num2str(round(time21(end),2),'%.1f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\returnFuel', namelist{j},returnparam ,'}{ ', num2str(round(mFuel22(1),1),'%.1f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\seconddExergy', namelist{j},returnparam ,'}{ ', num2str(dExergy_21(end)/10^9,'%.2f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\secondSpecExergy', namelist{j},returnparam ,'}{ ', num2str(Spec_Exergy_21/10^6,'%.2f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\secondExergyEff', namelist{j},returnparam ,'}{ ', num2str(eff_Exergy_2(end)*100,'%.3f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\secondExergythirdStage', namelist{j},returnparam ,'}{ ', num2str(dExergy_21_3(end)/10^9,'%.3f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\secondEnergy', namelist{j},returnparam ,'}{ ', num2str((Exergy_tot21)/10^9,'%.2f') , '}'), '-append','delimiter','','newline', 'pc');

dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\secondExergyDestroyed', namelist{j},returnparam ,'}{ ', num2str((Exergy_tot2-dExergy_21(end))/10^9,'%.3f') , '}'), '-append','delimiter','','newline', 'pc');

dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\secondWnextStage', namelist{j},returnparam ,'}{ ', num2str(W_21nextstage_pc(end) ,'%.3f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\secondWpayload', namelist{j},returnparam ,'}{ ', num2str(W_21Payload_pc(end) ,'%.3f') , '}'), '-append','delimiter','','newline', 'pc');


dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\SecondDist', namelist{j},returnparam ,'}{ ', num2str(round(total_dist_cov21,1),'%.1f') , '}'), '-append','delimiter','','newline', 'pc');

if auxdata.returnMode ==1
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\returnEnergy', namelist{j},returnparam ,'}{ ', num2str((Exergy_tot22)/10^9,'%.2f') , '}'), '-append','delimiter','','newline', 'pc');
end
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\returnDist', namelist{j},returnparam ,'}{ ', num2str(round(total_dist_cov22,1),'%.1f') , '}'), '-append','delimiter','','newline', 'pc');

dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\WDsecond', namelist{j},returnparam ,'}{ ', num2str(W_D21_pc(end) ,'%.2f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\Plosssecond', namelist{j},returnparam ,'}{ ', num2str(P_loss21_pc(end) ,'%.2f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\WmFsecond', namelist{j},returnparam ,'}{ ', num2str(W_mF21_pc(end) ,'%.2f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\Wsecond', namelist{j},returnparam ,'}{ ', num2str(W_2_pc(end) ,'%.2f') , '}'), '-append','delimiter','','newline', 'pc');

dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\PlosssecondCombined', namelist{j},returnparam ,'}{ ', num2str(P_loss21_pc(end)+W_mF21_pc(end) ,'%.2f') , '}'), '-append','delimiter','','newline', 'pc');
if auxdata.returnMode ==1
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\WDreturn', namelist{j},returnparam ,'}{ ', num2str(W_D22_pc(end) ,'%.2f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\PlossreturnCombined', namelist{j},returnparam ,'}{ ', num2str(P_loss22_pc(end)+W_mF22_pc(end) ,'%.2f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\WmFreturn', namelist{j},returnparam ,'}{ ', num2str(W_mF22_pc(end) ,'%.2f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\Wreturn', namelist{j},returnparam ,'}{ ', num2str(W_22_pc(end) ,'%.2f') , '}'), '-append','delimiter','','newline', 'pc');
end

qlt20 = find(q3<5000);
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\thirdqOverFive', namelist{j},returnparam ,'}{ ', num2str(round(time3(qlt20(1))-time3(1),2),'%.1f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\thirdmaxAoA', namelist{j},returnparam ,'}{ ', num2str(rad2deg(max(aoa3)),'%.1f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\thirdcircv', namelist{j},returnparam ,'}{ ', num2str(v3exo(end),'%.1f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\thirdcircm', namelist{j},returnparam ,'}{ ', num2str(mexo(end),'%.1f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\thirdmFuel', namelist{j},returnparam ,'}{ ', num2str(auxdata.Stage3.mTot-m3_4,'%.1f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\thirddExergy', namelist{j},returnparam ,'}{ ', num2str(dExergy_3/10^9,'%.2f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\thirddExergyEff', namelist{j},returnparam ,'}{ ', num2str(eff_Exergy_3*100,'%.3f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\thirdEnergy', namelist{j},returnparam ,'}{ ', num2str((Exergy_tot3)/10^9,'%.2f') , '}'), '-append','delimiter','','newline', 'pc');


dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\WDthree', namelist{j},returnparam ,'}{ ', num2str(W_D3_pc(end) ,'%.2f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\Wthree', namelist{j},returnparam ,'}{ ', num2str(W_3_pc(end) ,'%.2f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\WmFthree', namelist{j},returnparam ,'}{ ', num2str(W_mF3_pc(end) ,'%.2f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\HTloss', namelist{j},returnparam ,'}{ ', num2str(HT_loss_pc(end) ,'%.2f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\WHSthree', namelist{j},returnparam ,'}{ ', num2str(W_HS3_pc(end) ,'%.2f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\Plossthree ', namelist{j},returnparam ,'}{ ', num2str(P_loss3_pc(end) ,'%.2f') , '}'), '-append','delimiter','','newline', 'pc');

dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\PlossthreeCombined ', namelist{j},returnparam ,'}{ ', num2str(P_loss3_pc(end)+ W_mF3_pc(end)+HT_loss_pc(end),'%.2f') , '}'), '-append','delimiter','','newline', 'pc');


dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\thirdExergyDestroyed', namelist{j},returnparam ,'}{ ', num2str((Exergy_tot3-dExergy_3(end))/10^9,'%.3f') , '}'), '-append','delimiter','','newline', 'pc');



dlmwrite(strcat('LatexInputs.txt'),strcat('\newcommand{\totaldExergy', namelist{j},returnparam ,'}{ ', num2str((dExergy_1(end)+dExergy_21(end)+dExergy_3)/10^9,'%.2f') , '}'), '-append','delimiter','','newline', 'pc');


% if auxdata.returnMode == 1
if length(output)>1
ResultsMatrix(:,j) = [ThirdStagePayloadMass; alt21(1)/1000; v21(1); rad2deg(gamma21(1));...
    alt21(end)/1000; v21(end); rad2deg(gamma21(end)); q21(end)/1000; LD21(end);...
    time21(end); time3(qlt20(1))-time3(1); max(aoa3); v3exo(end); m3(end); mFuel22(1);...
    dExergy_1(end)/10^9; dExergy_21(end)/10^9; dExergy_3/10^9; eff_Exergy_1(end)*100;...
    eff_Exergy_2(end)*100; eff_Exergy_3*100; auxdata.Stage3.mTot-m3_4;eff_Exergy_total; (ThirdStagePayloadMass + output{3}.result.objective)/-output{3}.result.objective*100; total_dist_cov21;total_dist_cov22; Exergy_tot1-dExergy_1(end); Exergy_tot21-dExergy_21(end); Exergy_tot3-dExergy_3(end)]; % for variation study
else
 ResultsMatrix(:,j) = [ThirdStagePayloadMass; alt21(1)/1000; v21(1); rad2deg(gamma21(1));...
    alt21(end)/1000; v21(end); rad2deg(gamma21(end)); q21(end)/1000; LD21(end);...
    time21(end); time3(qlt20(1))-time3(1); max(aoa3); v3exo(end); m3(end); mFuel22(1);...
    dExergy_1(end)/10^9; dExergy_21(end)/10^9; dExergy_3/10^9; eff_Exergy_1(end)*100;...
    eff_Exergy_2(end)*100; eff_Exergy_3*100; auxdata.Stage3.mTot-m3_4;eff_Exergy_total; 0; total_dist_cov21;total_dist_cov22; Exergy_tot1-dExergy_1(end); Exergy_tot21-dExergy_21(end); Exergy_tot3-dExergy_3(end)]; % for variation study
   
end


% else
%  ResultsMatrix(:,j) = [ThirdStagePayloadMass;alt21(1)/1000;v21(1);rad2deg(gamma21(1));alt21(end)/1000;v21(end);rad2deg(gamma21(end));q21(end)/1000;LD21(end);time21(end);time3(qlt20(1))-time3(1);max(aoa3);v3exo(end);m3(end);mFuel22(1)]; % for variation study
%    
% end
% if auxdata.returnMode == 1
%     ResultsMatrix(end+1,j) = mFuel22(1);
% end
% movefile(strcat('LatexInputs.txt'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));


%% Bound Check
% Peform check to see if any of the states are hitting their bounds. This
% is an error if the bound is not intended to constrain the state. Fuel
% mass and throttle are not checked, as these will always hit bounds. 

for i = 1: length(output{j}.result.solution.phase(2).state(1,:))
    if any(output{j}.result.solution.phase(1).state(:,i) == bounds.phase(1).state.lower(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 1 is hitting lower bound'))
    end
    
    if any(output{j}.result.solution.phase(1).state(:,i) == bounds.phase(1).state.upper(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 1 is hitting upper bound'))
    end
end

for i = 1: length(output{j}.result.solution.phase(2).state(1,:))-1
    if any(output{j}.result.solution.phase(2).state(:,i) == bounds.phase(2).state.lower(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 2 is hitting lower bound'))
    end
    
    if any(output{j}.result.solution.phase(2).state(:,i) == bounds.phase(2).state.upper(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 2 is hitting upper bound'))
    end
end

for i = 1: length(output{j}.result.solution.phase(3).state(1,:))-2
    if any(output{j}.result.solution.phase(3).state(:,i) == bounds.phase(3).state.lower(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 3 is hitting lower bound'))
    end
    
    if any(output{j}.result.solution.phase(3).state(:,i) == bounds.phase(3).state.upper(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 3 is hitting upper bound'))
    end
end

% Angle of attack is not checked on third stage, because angle of attack is hard constrained and should be checked manually. 
for i = [1:3 6: length(output{j}.result.solution.phase(4).state(1,:))]
    if any(output{j}.result.solution.phase(4).state(:,i) == bounds.phase(4).state.lower(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 4 is hitting lower bound'))
    end
    
    if any(output{j}.result.solution.phase(4).state(:,i) == bounds.phase(4).state.upper(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 4 is hitting upper bound'))
    end
end

%% Plot combined trajectory




        figure(2301)
           fig = gcf;
        set(fig,'Position',[200 0 1200 1000])
hold on
% axesm('pcarree','MapLatLimit',[min(rad2deg(lat22))-5 max(rad2deg(latexo))+5],'MapLonLimit',[min(rad2deg(lon22))-5 max(rad2deg(lon21))+5])

if auxdata.mode == 0
    
worldmap([min(rad2deg(latexo))-4 max(rad2deg(lat1))+2],[min(rad2deg(lon21))-4 max(rad2deg(lon21))+4])
else

if returnMode == 1
worldmap([min(rad2deg(lat22))-4 max(rad2deg(latexo))+2],[min(rad2deg(lon22))-4 max(rad2deg(lon21))+4])
else
worldmap([min(rad2deg(lat21))-4 max(rad2deg(latexo))+2],[min(rad2deg(lon21))-4 max(rad2deg(lon21))+4])

end
end
northarrow('latitude',max(rad2deg(latexo))+0.5,'longitude',max(rad2deg(lon21))+3,'scaleratio',0.05)

geoshow('landareas.shp','FaceColor',[0.8 .8 0.8])
tightmap
mlabel('south')
if auxdata.mode == 0
   view([28 35]) 
else
view([-28 35])
end


% plotm(rad2deg(lat),rad2deg(lon+lon0))
a = plot3m(rad2deg(lat1),rad2deg(lon1),alt1*5,'color',[0 0.8 0],'LineWidth',1.7);
b = plot3m(rad2deg(lat21),rad2deg(lon21),alt21*5,'r','LineWidth',1.7);
if returnMode == 1
c = plot3m(rad2deg(lat22),rad2deg(lon22),alt22*5,'color',[1 0 1],'LineWidth',1.7);
end
d = plot3m(rad2deg(lat3),rad2deg(lon3)',alt3*5,'b','LineWidth',1.7);
e = plot3m(rad2deg(latexo),rad2deg(lonexo),altexo*5,'c','LineWidth',1.7);



 ht = text(-500000,-1600000,'Australia');
 set(ht,'Rotation',20)

plotm(rad2deg(lat1),rad2deg(lon1),'color',[0 0.8 0])
plotm(rad2deg(lat21),rad2deg(lon21),'r')
if returnMode == 1
plotm(rad2deg(lat22),rad2deg(lon22),'color',[1 0 1])
end
plotm(rad2deg(lat3),rad2deg(lon3)','b')
plotm(rad2deg(latexo),rad2deg(lonexo),'c')
   
    cities = shaperead('worldcities', 'UseGeoCoords', true);
lats = extractfield(cities,'Lat');
lons = extractfield(cities,'Lon');
geoshow(lats, lons,...
        'DisplayType', 'point',...
        'Marker', 'o',...
        'MarkerEdgeColor', 'r',...
        'MarkerFaceColor', 'r',...
        'MarkerSize', 2)

    % Plot vertical lines
for i = 0:20/(time1(end)-time1(1)):2 % goes to 2 because taking away remainder causes time to go negative. extrapolations do nothing
    time1_temp = (time1(end)-time1(1))*i + time1(1);
    alt1_temp = interp1(time1,alt1,time1_temp);
    lat1_temp = interp1(time1,lat1,time1_temp);
    lon1_temp = interp1(time1,lon1,time1_temp);

    plot3m(rad2deg([lat1_temp  ,lat1_temp ]), rad2deg([lon1_temp,  lon1_temp]),[0,alt1_temp*5],'color',[0.3 0.3 0.3]);

end
for i = 0:20/(time21(end)-time21(1)):2
    time21_temp = (time21(end)-time21(1))*i + time21(1) - rem(time21(1),20);
    alt21_temp = interp1(time21,alt21,time21_temp);
    lat21_temp = interp1(time21,lat21,time21_temp);
    lon21_temp = interp1(time21,lon21,time21_temp);

    plot3m(rad2deg([lat21_temp  ,lat21_temp ]), rad2deg([lon21_temp,  lon21_temp]),[0,alt21_temp*5],'color',[0.3 0.3 0.3]);

end
if returnMode == 1
for i = 0:20/(time22(end)-time22(1)):2
    time22_temp = (time22(end)-time22(1))*i + time22(1) - rem(time22(1),20) - rem(time21(1),20);
    alt22_temp = interp1(time22,alt22,time22_temp);
    lat22_temp = interp1(time22,lat22,time22_temp);
    lon22_temp = interp1(time22,lon22,time22_temp);

    plot3m(rad2deg([lat22_temp  ,lat22_temp ]), rad2deg([lon22_temp,  lon22_temp]),[0,alt22_temp*5],'color',[0.3 0.3 0.3]);

end
end
for i = 0:20/(time3(end)-time3(1)):2
    time3_temp = (time3(end)-time3(1))*i + time3(1) - rem(time3(1),20) - rem(time22(1),20) - rem(time21(1),20);
    alt3_temp = interp1(time3,alt3,time3_temp);
    lat3_temp = interp1(time3,lat3,time3_temp);
    lon3_temp = interp1(time3,lon3,time3_temp);
    plot3m(rad2deg([lat3_temp  ,lat3_temp ]), rad2deg([lon3_temp,  lon3_temp]),[0,alt3_temp*5],'color',[0.3 0.3 0.3]);

end
for i = 0:20/(timeexo(end)-timeexo(1)):2
    timeexo_temp = (timeexo(end)-timeexo(1))*i + timeexo(1) - rem(timeexo(1),20) - rem(time3(1),20) - rem(time22(1),20) - rem(time21(1),20);
    altexo_temp = interp1(timeexo,altexo,timeexo_temp);
    latexo_temp = interp1(timeexo,latexo,timeexo_temp);
    lonexo_temp = interp1(timeexo,lonexo,timeexo_temp);
    plot3m(rad2deg([latexo_temp  ,latexo_temp ]), rad2deg([lonexo_temp,  lonexo_temp]),[0,altexo_temp*5],'color',[0.3 0.3 0.3]);
end

scaleruler on
setm(handlem('scaleruler1'), ...
    'XLoc',-9.5e4,'YLoc',-1.5e6,...
    'MajorTick',0:200:400,'TickDir','down','RulerStyle','patches')

if returnMode == 1
legend([a b c d e],{'First Stage Ascent','Second Stage Ascent', 'Second Stage Return','Third Stage Powered Ascent','Third Stage Unpowered Ascent'},'Location','east');
else
 legend([a b d e],{'First Stage Ascent','Second Stage Ascent','Third Stage Powered Ascent','Third Stage Unpowered Ascent'},'Location','east');
end
% zoom(30) 
%% Plot Visualisation of Net ISP & other performance
plotperform = 'no'

if strcmp(plotperform,'yes')
    if auxdata.returnMode == 0 && auxdata.mode == 1
addpath('C:\Users\uqsforb1\Documents\GitHub\LODESTAR_FINAL\Results\mode900')
    
ConstqStates = dlmread('SecondStageStatesConstq'); % import constant dynamic pressure rsult for comparison
ConstqStates(ConstqStates(:,12)==0,:) = []; % remove return data
    end
    
% MList = [5:0.1:10];
altlist = 24000:100:35000 ;
alphalist = 0:.1:5;

[alphagrid,altgrid] = ndgrid(alphalist,altlist);
figure(501)
title('Net Isp')
hold on

figure(502)
title('Thrust')
hold on
for temp = 1:5
    
% v_temp = 1620+300*(temp-1);
% M_temp = 6.5+temp/4
M_temp = 4.5+temp
for i = 1:numel(alphagrid)
        I = cell(1, ndims(alphagrid)); 
    [I{:}] = ind2sub(size(alphagrid),i);
    alpha_temp = alphagrid(I{1},I{2});
    alt_temp = altgrid(I{1},I{2});
    
    v_temp = M_temp*300; % velocity is only approximate
    
[~,~,~,~,~,~, ~, ~, Fdgrid(I{(1)},I{(2)}), ~,Lgrid(I{(1)},I{(2)}),Fueldtgrid(I{(1)},I{(2)}),Tgrid(I{(1)},I{(2)}),Ispgrid(I{(1)},I{(2)}),~,~,~,~,T1grid(I{(1)},I{(2)}),P1grid(I{(1)},I{(2)}),M1grid(I{(1)},I{(2)})] = SPARTANDynamics(0, alt_temp, v_temp,auxdata,0,0,0,deg2rad(alpha_temp),0,1, mFuel21(1),mFuel21(1),mFuel21(end), 1, 0);

end

% [~,~,~,~,~,~, ~, ~, Fd, ~,L,Fueldt,T,Isp,~,~,~,~] = SPARTANDynamics(0, 34000, 1600,auxdata,0,0,0,deg2rad(1),0,1, mFuel21(1),mFuel21(1),mFuel21(end), 1, 0)

figure(501)
subplot(2,3,temp);
hold on
colormap jet
contourf(alphagrid,altgrid/1000,(Tgrid-Fdgrid)./Fueldtgrid/9.81,500,'LineColor','none')
xlabel('Angle of Attack (deg)','FontSize',8);
ylabel('Altitude (km)','FontSize',8);
title(['Mach ' num2str(M_temp)],'FontSize',10)
set(gca,'XTick',0:5)
caxis([0 1500]);
alt_temp = interp1(M21,alt21,M_temp);
alpha_temp = interp1(M21,alpha21,M_temp);
plot(rad2deg(alpha_temp),alt_temp/1000,'o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','m')
if auxdata.returnMode == 0 && auxdata.mode == 1
alt_temp_constq = interp1(ConstqStates(:,3),ConstqStates(:,2),M_temp);
alpha_temp_constq = interp1(ConstqStates(:,3),ConstqStates(:,9),M_temp);
plot(rad2deg(alpha_temp_constq ),alt_temp_constq /1000,'o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','c')

legend('-','Opimised Trajectory','Constant q Trajectory')
end

colorbar
h = colorbar;
ylabel(h, 'Net Isp')

mean_isp68 = mean(interp1(M21,Isp21,6:0.01:8.3))

mean_netisp68 = mean(interp1(M21,IspNet21,6:0.01:8.3))

time68 = interp1(M21,time21,8)-interp1(M21,time21,6)

figure(503)
subplot(2,3,temp);
hold on
colormap jet
contourf(alphagrid,altgrid/1000,Tgrid,500,'LineColor','none')
xlabel('Angle of Attack (deg)','FontSize',8);
ylabel('Altitude (km)','FontSize',8);
title(['Mach ' num2str(M_temp)],'FontSize',10)
set(gca,'XTick',0:5)
caxis([0 150000]);
alt_temp = interp1(M21,alt21,M_temp);
alpha_temp = interp1(M21,alpha21,M_temp);
plot(rad2deg(alpha_temp),alt_temp/1000,'o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','m')


figure(504)
subplot(2,3,temp);
hold on
colormap jet
contourf(alphagrid,altgrid/1000,Fueldtgrid,500,'LineColor','none')
xlabel('Angle of Attack (deg)','FontSize',8);
ylabel('Altitude (km)','FontSize',8);
title(['Mach ' num2str(M_temp)],'FontSize',10)
set(gca,'XTick',0:5)
% caxis([0 150000]);
alt_temp = interp1(M21,alt21,M_temp);
alpha_temp = interp1(M21,alpha21,M_temp);
plot(rad2deg(alpha_temp),alt_temp/1000,'o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','m')

figure(505)
subplot(3,3,temp);
hold on
colormap jet
contourf(alphagrid,altgrid/1000,Fdgrid,500,'LineColor','none')
xlabel('Angle of Attack (deg)','FontSize',8);
ylabel('Altitude (km)','FontSize',8);
title(['Mach ' num2str(M_temp)],'FontSize',10)
set(gca,'XTick',0:5)
% caxis([0 150000]);
alt_temp = interp1(M21,alt21,M_temp);
alpha_temp = interp1(M21,alpha21,M_temp);
plot(rad2deg(alpha_temp),alt_temp/1000,'o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','m')

end


saveas(figure(501),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('NetIsp',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(501),strcat('NetIsp',namelist{j}),'-dpng');
movefile(strcat('NetIsp',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));
saveas(figure(503),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('Thrust',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(503),strcat('Thrust',namelist{j}),'-dpng');
movefile(strcat('Thrust',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));
saveas(figure(504),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('Fueldt',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(504),strcat('Fueldt',namelist{j}),'-dpng');
movefile(strcat('Fueldt',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));

close(figure(501))
close(figure(503))
close(figure(504))


saveas(figure(505),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('Fd',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(505),strcat('Fd',namelist{j}),'-dpng');
movefile(strcat('Fd',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));

close(figure(505))

for temp = 1:6
    
% v_temp = 1620+300*(temp-1);
M_temp = 6.5+temp/4

for i = 1:numel(alphagrid)
        I = cell(1, ndims(alphagrid)); 
    [I{:}] = ind2sub(size(alphagrid),i);
    alpha_temp = alphagrid(I{1},I{2});
    alt_temp = altgrid(I{1},I{2});
    
    v_temp = M_temp*300; % velocity is only approximate
    
[~,~,~,~,~,~, ~, ~, Fdgrid(I{(1)},I{(2)}), ~,Lgrid(I{(1)},I{(2)}),Fueldtgrid(I{(1)},I{(2)}),Tgrid(I{(1)},I{(2)}),Ispgrid(I{(1)},I{(2)}),~,~,~,~,T1grid(I{(1)},I{(2)}),P1grid(I{(1)},I{(2)}),M1grid(I{(1)},I{(2)})] = SPARTANDynamics(0, alt_temp, v_temp,auxdata,0,0,0,deg2rad(alpha_temp),0,1, mFuel21(1),mFuel21(1),mFuel21(end), 1, 0);

end


figure(506)
subplot(2,3,temp);
hold on
colormap jet
contourf(alphagrid,altgrid/1000,Ispgrid,500,'LineColor','none')
xlabel('Angle of Attack (deg)','FontSize',8);
ylabel('Altitude (km)','FontSize',8);
title(['Isp-Mach ' num2str(M_temp)],'FontSize',10)
set(gca,'XTick',0:5)
% caxis([0 1500]);
alt_temp = interp1(M21,alt21,M_temp);
alpha_temp = interp1(M21,alpha21,M_temp);
plot(rad2deg(alpha_temp),alt_temp/1000,'o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','m')
if auxdata.returnMode == 0 && auxdata.mode == 1
alt_temp_constq = interp1(ConstqStates(:,3),ConstqStates(:,2),M_temp);
alpha_temp_constq = interp1(ConstqStates(:,3),ConstqStates(:,9),M_temp);
plot(rad2deg(alpha_temp_constq ),alt_temp_constq /1000,'o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','c')

legend('-','Opimised Trajectory','Constant q Trajectory')
end

colorbar
h = colorbar;
ylabel(h, 'Isp')

figure(507)
subplot(2,3,temp);
hold on
colormap jet
contourf(alphagrid,altgrid/1000,T1grid,500,'LineColor','none')
xlabel('Angle of Attack (deg)','FontSize',8);
ylabel('Altitude (km)','FontSize',8);
title(['T1 M ' num2str(M_temp)],'FontSize',10)
set(gca,'XTick',0:5)
% caxis([0 1500]);
alt_temp = interp1(M21,alt21,M_temp);
alpha_temp = interp1(M21,alpha21,M_temp);
plot(rad2deg(alpha_temp),alt_temp/1000,'o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','m')


figure(508)
subplot(2,3,temp);
hold on
colormap jet
contourf(alphagrid,altgrid/1000,P1grid,500,'LineColor','none')
xlabel('Angle of Attack (deg)','FontSize',8);
ylabel('Altitude (km)','FontSize',8);
title(['P1 M ' num2str(M_temp)],'FontSize',10)
set(gca,'XTick',0:5)
% caxis([0 1500]);
alt_temp = interp1(M21,alt21,M_temp);
alpha_temp = interp1(M21,alpha21,M_temp);
plot(rad2deg(alpha_temp),alt_temp/1000,'o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','m')

figure(509)
subplot(2,3,temp);
hold on
colormap jet
contourf(alphagrid,altgrid/1000,M1grid,500,'LineColor','none')
xlabel('Angle of Attack (deg)','FontSize',8);
ylabel('Altitude (km)','FontSize',8);
title(['M1 M ' num2str(M_temp)],'FontSize',10)
set(gca,'XTick',0:5)
% caxis([0 1500]);
alt_temp = interp1(M21,alt21,M_temp);
alpha_temp = interp1(M21,alpha21,M_temp);
plot(rad2deg(alpha_temp),alt_temp/1000,'o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','m')

end



saveas(figure(506),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('Isp-Mach',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(506),strcat('Isp-Mach',namelist{j}),'-dpng');
movefile(strcat('Isp-Mach',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));



saveas(figure(507),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('T1',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(507),strcat('T1',namelist{j}),'-dpng');
movefile(strcat('T1',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));


close(figure(506))
close(figure(507))

saveas(figure(508),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('P1',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(508),strcat('P1',namelist{j}),'-dpng');
movefile(strcat('P1',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));
saveas(figure(509),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('M1',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(509),strcat('M1',namelist{j}),'-dpng');
movefile(strcat('M1',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));

close(figure(508))
close(figure(509))

%% Plot End of ascent net Isp comparison
% 
% 
% % MList = [5:0.1:10];
% altlist = 35000:100:45000 ;
% alphalist = 5:.1:10;
% 
% [alphagrid,altgrid] = ndgrid(alphalist,altlist);
% figure(502)
% title('Net Isp')
% hold on
% for temp = 1:5
%     
% v_temp = interp1(time21,v21,time21(end)-21+5*(temp-1));
% 
% for i = 1:numel(alphagrid)
%         I = cell(1, ndims(alphagrid)); 
%     [I{:}] = ind2sub(size(alphagrid),i);
%     alpha_temp = alphagrid(I{1},I{2});
%     alt_temp = altgrid(I{1},I{2});
%     
%     
% [~,~,~,~,~,~, ~, ~, Fdgrid(I{(1)},I{(2)}), ~,Lgrid(I{(1)},I{(2)}),Fueldtgrid(I{(1)},I{(2)}),Tgrid(I{(1)},I{(2)}),Ispgrid(I{(1)},I{(2)}),~,~,~,~] = SPARTANDynamics(0, alt_temp, v_temp,auxdata,0,0,0,deg2rad(alpha_temp),0,1, mFuel21(1),mFuel21(1),mFuel21(end), 1, 0);
% 
% end
% 
% % [~,~,~,~,~,~, ~, ~, Fd, ~,L,Fueldt,T,Isp,~,~,~,~] = SPARTANDynamics(0, 34000, 1600,auxdata,0,0,0,deg2rad(1),0,1, mFuel21(1),mFuel21(1),mFuel21(end), 1, 0)
% 
% 
% subplot(2,3,temp);
% hold on
% colormap jet
% contourf(alphagrid,altgrid/1000,(Tgrid-Fdgrid)./Fueldtgrid/9.81,500,'LineColor','none')
% xlabel('Angle of Attack (deg)','FontSize',8);
% ylabel('Altitude (km)','FontSize',8);
% title([num2str(v_temp) ' m/s'],'FontSize',10)
% set(gca,'XTick',5:10)
% caxis([-1000 500]);
% alt_temp = interp1(time21,alt21,time21(end)-21+5*(temp-1));
% alpha_temp = interp1(time21,alpha21,time21(end)-21+5*(temp-1));
% plot(rad2deg(alpha_temp),alt_temp/1000,'o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','m')
% 
% end

%%

c = colorbar;
set(c, 'Position', [.70 .13 .0581 .3])
ylabel(c, 'Net Isp')

% figure(503)
% contourf(alphagrid,altgrid,Fdgrid,1000,'LineWidth',0.)
% c = colorbar;
% 
% figure(504)
% contourf(alphagrid,altgrid,Tgrid,1000,'LineWidth',0.)
% title('Thrust')
% c = colorbar;
% 
% figure(505)
% contourf(alphagrid,altgrid,Fueldtgrid,1000,'LineWidth',0.)
% c = colorbar;



end

%% SAVE FIGS

saveas(figure(311),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('ThirdStage',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(311),strcat('ThirdStage',namelist{j}),'-dpng');
movefile(strcat('ThirdStage',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));
saveas(figure(3223),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('AoAComparison',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(3223),strcat('AoAComparison',namelist{j}),'-dpng');
movefile(strcat('AoAComparison',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));
saveas(figure(211),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('SecondStage',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(211),strcat('SecondStage',namelist{j}),'-dpng');
movefile(strcat('SecondStage',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));
saveas(figure(221),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('Return',namelist{j},'.fig')]); 
set(gcf, 'PaperPositionMode', 'auto');
print(figure(221),strcat('Return',namelist{j}),'-dpng');
movefile(strcat('Return',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));
saveas(figure(701),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('HeatFlux',namelist{j},'.fig')]); 
set(gcf, 'PaperPositionMode', 'auto');
print(figure(701),strcat('HeatFlux',namelist{j}),'-dpng');
movefile(strcat('HeatFlux',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));
saveas(figure(2410),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('Hamiltonian',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(2410),strcat('Hamiltonian',namelist{j}),'-dpng');
movefile(strcat('Hamiltonian',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));
saveas(figure(2420),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('Verification',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(2420),strcat('Verification',namelist{j}),'-dpng');
movefile(strcat('Verification',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));
saveas(figure(550),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('ForwardError',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(550),strcat('ForwardError',namelist{j}),'-dpng');
movefile(strcat('ForwardError',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));
saveas(figure(212),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('Forward1',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(212),strcat('Forward1',namelist{j}),'-dpng');
movefile(strcat('Forward1',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));
saveas(figure(213),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('Forward2',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(213),strcat('Forward2',namelist{j}),'-dpng');
movefile(strcat('Forward2',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));
% saveas(figure(2100),[sprintf('../ArchivedResults/%s',Timestamp),filesep,'eq.fig']);
% saveas(figure(2110),[sprintf('../ArchivedResults/%s',Timestamp),filesep,'ISP.fig']);
if IspPlot == 1
saveas(figure(2100),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('eq',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(2100),strcat('eq',namelist{j}),'-dpng');
movefile(strcat('eq',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));
saveas(figure(2110),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('Isp',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(2110),strcat('Isp',namelist{j}),'-dpng');
movefile(strcat('Isp',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));
saveas(figure(2210),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('returnIsp',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(2210),strcat('returnIsp',namelist{j}),'-dpng');
movefile(strcat('returnIsp',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));
end
saveas(figure(2301),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('GroundTrack',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperPositionMode', 'auto');
print(figure(2301),strcat('GroundTrack',namelist{j}),'-dpng');
movefile(strcat('GroundTrack',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));
% saveas(figure(501),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('NetIsp',namelist{j},'.fig')]);
% set(gcf, 'PaperPositionMode', 'auto');
% print(figure(501),strcat('NetIsp',namelist{j}),'-dpng');
% movefile(strcat('NetIsp',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));
% saveas(figure(503),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('Thrust',namelist{j},'.fig')]);
% set(gcf, 'PaperPositionMode', 'auto');
% print(figure(503),strcat('Thrust',namelist{j}),'-dpng');
% movefile(strcat('Thrust',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));
% saveas(figure(504),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('Fueldt',namelist{j},'.fig')]);
% set(gcf, 'PaperPositionMode', 'auto');
% print(figure(504),strcat('Fueldt',namelist{j}),'-dpng');
% movefile(strcat('Fueldt',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));
% saveas(figure(505),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('Fd',namelist{j},'.fig')]);
% set(gcf, 'PaperPositionMode', 'auto');
% print(figure(505),strcat('Fd',namelist{j}),'-dpng');
% movefile(strcat('Fd',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));
% saveas(figure(506),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('Isp-Mach',namelist{j},'.fig')]);
% set(gcf, 'PaperPositionMode', 'auto');
% print(figure(506),strcat('Isp-Mach',namelist{j}),'-dpng');
% movefile(strcat('Isp-Mach',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));
% 
% saveas(figure(507),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('T1',namelist{j},'.fig')]);
% set(gcf, 'PaperPositionMode', 'auto');
% print(figure(507),strcat('T1',namelist{j}),'-dpng');
% movefile(strcat('T1',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));
% saveas(figure(508),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('P1',namelist{j},'.fig')]);
% set(gcf, 'PaperPositionMode', 'auto');
% print(figure(508),strcat('P1',namelist{j}),'-dpng');
% movefile(strcat('P1',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));
% saveas(figure(509),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('M1',namelist{j},'.fig')]);
% set(gcf, 'PaperPositionMode', 'auto');
% print(figure(509),strcat('M1',namelist{j}),'-dpng');
% movefile(strcat('M1',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));


saveas(figure(502),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('NetIspPullup',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(502),strcat('NetIspPullup',namelist{j}),'-dpng');
movefile(strcat('NetIspPullup',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));


close all


%% Overpressure
% calculated from Carlson and Meglieri, Review of Sonic-Boom Generation
% Theory and Prediction Methods
% https://asa.scitation.org/doi/pdf/10.1121/1.1912901
% https://asa.scitation.org/doi/pdf/10.1121/1.1912909

beta = (([M21; M22(M22>1)]).^2-1).^(0.5);
W = ([auxdata.Stage2.mStruct+auxdata.Stage3.mTot+mFuel21; auxdata.Stage2.mStruct+mFuel22(M22>1)])*2.20462;% Total Weight in lbs
gamma_op = 1.4;% Ratio of specific heats
P_h = [P021; P022(M22>1)]*20.885; % Pressure at altitude in psf 
l  = 22.94*3.28; % Reference length
h  = [alt21; alt22(M22>1)]*3.28; % altitude in ft

KA = h/120e3*0.2+1; % an estimated linear regression (will be a little off, but will make at max 5% difference)
P_g  = 2116.23; % Ground Pressure psf
KL = beta.*W./(gamma_op.*P_h.*[M21; M22(M22>1)].^2.*l.^2); % lift parameter
Ks = KL/0.03.*(0.135-0.045)+0.045; %interpolate for Ks (fig 12 in paper)


OP = 1.9.*beta.^(0.25) .*Ks.*KA.*sqrt(P_h.*P_g)./(h./l).^(0.75); %overpressure

figure(1022)
plot([time21;time22(M22>1)],OP*47.880172);

xlabel('Time (s)');
ylabel('Overpressure (Pa)')
title('Sonic Boom Ground Effects')

saveas(figure(1022),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,strcat('OverPressure',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(1022),strcat('OverPressure',namelist{j}),'-dpng');
movefile(strcat('OverPressure',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));

end
end

%% Create latex table



if auxdata.mode == 10
RegressionList = [-5 -2.5 0 2.5 5];
elseif auxdata.mode == 2 || auxdata.mode == 7
RegressionList = [-20 -10 0 10 20];  
elseif auxdata.mode == 5
RegressionList = [-80 -50 0 7 15]    ;
else
RegressionList = [-10 -5 0 5 10];
end

if mode == 1 || mode == 90
dlmwrite(strcat('LatexInputs.txt'),'\begin{tabular}{l c } ', '-append' , 'delimiter','','newline', 'pc')
else
 dlmwrite(strcat('LatexInputs.txt'),'\begin{tabular}{l c c c c c c} ', '-append' , 'delimiter','','newline', 'pc')
end

dlmwrite(strcat('LatexInputs.txt'),'\hline \textbf{Trajectory Condition}', '-append' , 'delimiter','','newline', 'pc')
for i = 1:length(namelist)
dlmwrite(strcat('LatexInputs.txt'), strcat('& ',namelist{i}) , '-append' , 'delimiter','','newline', 'pc')
end
dlmwrite(strcat('LatexInputs.txt'), '& $\Delta/\Delta$\%', '-append' , 'delimiter','','newline', 'pc')
dlmwrite(strcat('LatexInputs.txt'),'\\', '-append' , 'delimiter','','newline', 'pc')
%%
dlmwrite(strcat('LatexInputs.txt'),'\hline \textbf{Payload to Orbit (kg)}', '-append' , 'delimiter','','newline', 'pc')
for i = 1:length(namelist)
dlmwrite(strcat('LatexInputs.txt'),strcat('& \textbf{\PayloadToOrbit', namelist{i},returnparam,'}' ), '-append' , 'delimiter','','newline', 'pc')
end
if mode ~= 1 && mode ~= 90 && mode ~=0 && mode ~= 101 && mode ~=1000
    if issorted(ResultsMatrix(1,:)) || issorted(flip(ResultsMatrix(1,:)))
        Regression = RegressionList'\ResultsMatrix(1,:)';
        dlmwrite(strcat('LatexInputs.txt'),strcat('&\textbf{',num2str(round(Regression,1)),'}'), '-append' , 'delimiter','','newline', 'pc')
    else
        dlmwrite(strcat('LatexInputs.txt'),'& -', '-append' , 'delimiter','','newline', 'pc')
    end
end
dlmwrite(strcat('LatexInputs.txt'),'\\', '-append' , 'delimiter','','newline', 'pc')
%%
if length(output)>1
dlmwrite(strcat('LatexInputs.txt'),' \textbf{Payload Variation (\%)}', '-append' , 'delimiter','','newline', 'pc')
for i = 1:length(namelist)
dlmwrite(strcat('LatexInputs.txt'),strcat('& \PayloadVar', namelist{i},returnparam ), '-append' , 'delimiter','','newline', 'pc')
end
if mode ~= 1 && mode ~= 90 && mode ~=0 && mode ~= 101 && mode ~=1000
    if issorted(ResultsMatrix(24,:)) || issorted(flip(ResultsMatrix(24,:)))
        Regression = RegressionList'\ResultsMatrix(24,:)';
        dlmwrite(strcat('LatexInputs.txt'),strcat('&',num2str(round(Regression,2))), '-append' , 'delimiter','','newline', 'pc')
    else
        dlmwrite(strcat('LatexInputs.txt'),'& -', '-append' , 'delimiter','','newline', 'pc')
    end
end
dlmwrite(strcat('LatexInputs.txt'),'\\', '-append' , 'delimiter','','newline', 'pc')

end
%%
dlmwrite(strcat('LatexInputs.txt'),' \textbf{Total $\eta_{exergy}$ (\%)}', '-append' , 'delimiter','','newline', 'pc')
for i = 1:length(namelist)
dlmwrite(strcat('LatexInputs.txt'),strcat('& \textbf{\totalExergyEff', namelist{i},returnparam,'}' ), '-append','delimiter','','newline', 'pc')
end
if mode ~= 1 && mode ~= 90 && mode ~=0 && mode ~= 101 && mode ~=1000
    if issorted(ResultsMatrix(23,:)) || issorted(flip(ResultsMatrix(23,:)))
        Regression = RegressionList'\ResultsMatrix(23,:)';
        dlmwrite(strcat('LatexInputs.txt'),strcat('& \textbf{',num2str(round(Regression,5)),'}'), '-append' , 'delimiter','','newline', 'pc')
    else
        dlmwrite(strcat('LatexInputs.txt'),'& -', '-append' , 'delimiter','','newline', 'pc')
    end
end
dlmwrite(strcat('LatexInputs.txt'),'\\', '-append' , 'delimiter','','newline', 'pc')	
dlmwrite(strcat('LatexInputs.txt'),'\hline ', '-append' , 'delimiter','','newline', 'pc')
% %%
dlmwrite(strcat('LatexInputs.txt'),' \textbf{1$^{st}$ Stage $\eta_{exergy}$ (\%)}', '-append' , 'delimiter','','newline', 'pc')
for i = 1:length(namelist)
dlmwrite(strcat('LatexInputs.txt'),strcat('& \textbf{\firstExergyEff', namelist{i},returnparam,'}' ), '-append','delimiter','','newline', 'pc')
end
if mode ~= 1 && mode ~= 90 && mode ~=0 && mode ~= 101 && mode ~=1000
    if issorted(ResultsMatrix(19,:)) || issorted(flip(ResultsMatrix(19,:)))
        Regression = RegressionList'\ResultsMatrix(19,:)';
        dlmwrite(strcat('LatexInputs.txt'),strcat('& \textbf{',num2str(round(Regression,3)),'}'), '-append' , 'delimiter','','newline', 'pc')
    else
        dlmwrite(strcat('LatexInputs.txt'),'& -', '-append' , 'delimiter','','newline', 'pc')
    end
end
dlmwrite(strcat('LatexInputs.txt'),'\\', '-append' , 'delimiter','','newline', 'pc')
% %%
% dlmwrite(strcat('LatexInputs.txt'),' \textbf{1$^{st}$ Stage Destroyed Exergy (GJ)}', '-append' , 'delimiter','','newline', 'pc')
% for i = 1:length(namelist)
% dlmwrite(strcat('LatexInputs.txt'),strcat('& \textbf{\firstExergyDestroyed', namelist{i},returnparam,'}' ), '-append','delimiter','','newline', 'pc')
% end
% if mode ~= 1 && mode ~= 90 && mode ~=0
%     if issorted(ResultsMatrix(27,:)) || issorted(flip(ResultsMatrix(27,:)))
%         Regression = RegressionList'\ResultsMatrix(27,:)';
%         dlmwrite(strcat('LatexInputs.txt'),strcat('& \textbf{',num2str(round(Regression/10^9,3)),'}'), '-append' , 'delimiter','','newline', 'pc')
%     else
%         dlmwrite(strcat('LatexInputs.txt'),'& -', '-append' , 'delimiter','','newline', 'pc')
%     end
% end
% dlmwrite(strcat('LatexInputs.txt'),'\\', '-append' , 'delimiter','','newline', 'pc')

%%
dlmwrite(strcat('LatexInputs.txt'),' \textbf{Separation Alt, 1$\rightarrow$2 (km)}', '-append' , 'delimiter','','newline', 'pc')
for i = 1:length(namelist)
dlmwrite(strcat('LatexInputs.txt'),strcat('& \firstsecondSeparationAlt', namelist{i},returnparam ), '-append','delimiter','','newline', 'pc')
end
if mode ~= 1 && mode ~= 90 && mode ~=0 && mode ~= 101 && mode ~=1000
    if issorted(ResultsMatrix(2,:)) || issorted(flip(ResultsMatrix(2,:)))
        Regression = RegressionList'\ResultsMatrix(2,:)';
        dlmwrite(strcat('LatexInputs.txt'),strcat('&',num2str(round(Regression,2))), '-append' , 'delimiter','','newline', 'pc')
    else
        dlmwrite(strcat('LatexInputs.txt'),'& -', '-append' , 'delimiter','','newline', 'pc')
    end
end
dlmwrite(strcat('LatexInputs.txt'),'\\', '-append' , 'delimiter','','newline', 'pc')
%%
dlmwrite(strcat('LatexInputs.txt'),' \textbf{Separation v, 1$\rightarrow$2 (m/s)}', '-append' , 'delimiter','','newline', 'pc')
for i = 1:length(namelist)
dlmwrite(strcat('LatexInputs.txt'),strcat('& \firstsecondSeparationv', namelist{i} ,returnparam), '-append','delimiter','','newline', 'pc')
end
if mode ~= 1 && mode ~= 90 && mode ~=0 && mode ~= 101 && mode ~=1000
    if issorted(ResultsMatrix(3,:)) || issorted(flip(ResultsMatrix(3,:)))
        Regression = RegressionList'\ResultsMatrix(3,:)';
        dlmwrite(strcat('LatexInputs.txt'),strcat('&',num2str(round(Regression,2))), '-append' , 'delimiter','','newline', 'pc')
    else
        dlmwrite(strcat('LatexInputs.txt'),'& -', '-append' , 'delimiter','','newline', 'pc')
    end
end
dlmwrite(strcat('LatexInputs.txt'),'\\', '-append' , 'delimiter','','newline', 'pc')
%%
dlmwrite(strcat('LatexInputs.txt'),' \textbf{Separation $\gamma$, 1$\rightarrow$2 (deg)}', '-append' , 'delimiter','','newline', 'pc')
for i = 1:length(namelist)
dlmwrite(strcat('LatexInputs.txt'),strcat('& \firstsecondSeparationgamma', namelist{i},returnparam ), '-append','delimiter','','newline', 'pc')
end
if mode ~= 1 && mode ~= 90 && mode ~=0 && mode ~= 101 && mode ~=1000
    if issorted(ResultsMatrix(4,:)) || issorted(flip(ResultsMatrix(4,:)))
        Regression = RegressionList'\ResultsMatrix(4,:)';
        dlmwrite(strcat('LatexInputs.txt'),strcat('&',num2str(round(Regression,2))), '-append' , 'delimiter','','newline', 'pc')
    else
        dlmwrite(strcat('LatexInputs.txt'),'& -', '-append' , 'delimiter','','newline', 'pc')
    end
end
dlmwrite(strcat('LatexInputs.txt'),'\\', '-append' , 'delimiter','','newline', 'pc')
dlmwrite(strcat('LatexInputs.txt'),'\hline ', '-append' , 'delimiter','','newline', 'pc')
% 	%%
dlmwrite(strcat('LatexInputs.txt'),' \textbf{2$^{nd}$ Stage $\eta_{exergy}$ (\%)}', '-append' , 'delimiter','','newline', 'pc')
for i = 1:length(namelist)
dlmwrite(strcat('LatexInputs.txt'),strcat('& \textbf{\secondExergyEff', namelist{i},returnparam,'}' ), '-append','delimiter','','newline', 'pc')
end
if mode ~= 1 && mode ~= 90 && mode ~=0 && mode ~= 101 && mode ~=1000

    if issorted(ResultsMatrix(20,:)) || issorted(flip(ResultsMatrix(20,:)))
        Regression = RegressionList'\ResultsMatrix(20,:)';
        dlmwrite(strcat('LatexInputs.txt'),strcat('& \textbf{',num2str(round(Regression,3)),'}'), '-append' , 'delimiter','','newline', 'pc')
    else
        dlmwrite(strcat('LatexInputs.txt'),'& -', '-append' , 'delimiter','','newline', 'pc')
    end
end
dlmwrite(strcat('LatexInputs.txt'),'\\', '-append' , 'delimiter','','newline', 'pc')
% %%
% dlmwrite(strcat('LatexInputs.txt'),' \textbf{2$^{nd}$ Stage Destroyed Exergy (GJ)}', '-append' , 'delimiter','','newline', 'pc')
% for i = 1:length(namelist)
% dlmwrite(strcat('LatexInputs.txt'),strcat('& \textbf{\secondExergyDestroyed', namelist{i},returnparam,'}' ), '-append','delimiter','','newline', 'pc')
% end
% if mode ~= 1 && mode ~= 90 && mode ~=0
%     if issorted(ResultsMatrix(28,:)) || issorted(flip(ResultsMatrix(28,:)))
%         Regression = RegressionList'\ResultsMatrix(28,:)';
%         dlmwrite(strcat('LatexInputs.txt'),strcat('& \textbf{',num2str(round(Regression/10^9,3)),'}'), '-append' , 'delimiter','','newline', 'pc')
%     else
%         dlmwrite(strcat('LatexInputs.txt'),'& -', '-append' , 'delimiter','','newline', 'pc')
%     end
% end
% dlmwrite(strcat('LatexInputs.txt'),'\\', '-append' , 'delimiter','','newline', 'pc')

%%
dlmwrite(strcat('LatexInputs.txt'),' \textbf{Separation Alt, 2$\rightarrow$3 (km)}', '-append' , 'delimiter','','newline', 'pc')
for i = 1:length(namelist)
dlmwrite(strcat('LatexInputs.txt'),strcat('& \secondthirdSeparationAlt', namelist{i},returnparam), '-append','delimiter','','newline', 'pc');
end
if mode ~= 1 && mode ~= 90 && mode ~=0 && mode ~= 101 && mode ~=1000
    if issorted(ResultsMatrix(5,:)) || issorted(flip(ResultsMatrix(5,:)))
        Regression = RegressionList'\ResultsMatrix(5,:)';
        dlmwrite(strcat('LatexInputs.txt'),strcat('&',num2str(round(Regression,2))), '-append' , 'delimiter','','newline', 'pc')
    else
        dlmwrite(strcat('LatexInputs.txt'),'& -', '-append' , 'delimiter','','newline', 'pc')
    end
end
dlmwrite(strcat('LatexInputs.txt'),'\\', '-append' , 'delimiter','','newline', 'pc')
%%
dlmwrite(strcat('LatexInputs.txt'),' \textbf{Separation $v$, 2$\rightarrow$3 (m/s)}', '-append' , 'delimiter','','newline', 'pc')
for i = 1:length(namelist)
dlmwrite(strcat('LatexInputs.txt'),strcat('& \secondthirdSeparationv', namelist{i},returnparam ), '-append','delimiter','','newline', 'pc');
end
if mode ~= 1 && mode ~= 90 && mode ~=0 && mode ~= 101 && mode ~=1000
    if issorted(ResultsMatrix(6,:)) || issorted(flip(ResultsMatrix(6,:)))
        Regression = RegressionList'\ResultsMatrix(6,:)';
        dlmwrite(strcat('LatexInputs.txt'),strcat('&',num2str(round(Regression,2))), '-append' , 'delimiter','','newline', 'pc')
    else
        dlmwrite(strcat('LatexInputs.txt'),'& -', '-append' , 'delimiter','','newline', 'pc')
    end
end
dlmwrite(strcat('LatexInputs.txt'),'\\', '-append' , 'delimiter','','newline', 'pc')
%%
dlmwrite(strcat('LatexInputs.txt'),' \textbf{Separation $\gamma$, 2$\rightarrow$3 (deg)}', '-append' , 'delimiter','','newline', 'pc')
for i = 1:length(namelist)
dlmwrite(strcat('LatexInputs.txt'),strcat('& \secondthirdSeparationgamma', namelist{i},returnparam), '-append','delimiter','','newline', 'pc');
end
if mode ~= 1 && mode ~= 90 && mode ~=0 && mode ~= 101 && mode ~=1000
    if issorted(ResultsMatrix(7,:)) || issorted(flip(ResultsMatrix(7,:)))
        Regression = RegressionList'\ResultsMatrix(7,:)';
        dlmwrite(strcat('LatexInputs.txt'),strcat('&',num2str(round(Regression,2))), '-append' , 'delimiter','','newline', 'pc')
    else
        dlmwrite(strcat('LatexInputs.txt'),'& -', '-append' , 'delimiter','','newline', 'pc')
    end
end
dlmwrite(strcat('LatexInputs.txt'),'\\', '-append' , 'delimiter','','newline', 'pc')

%%
dlmwrite(strcat('LatexInputs.txt'),' \textbf{2$^{nd}$ Stage Flight Time (s)}', '-append' , 'delimiter','','newline', 'pc')
for i = 1:length(namelist)
dlmwrite(strcat('LatexInputs.txt'),strcat('& \secondFlightTime', namelist{i},returnparam), '-append','delimiter','','newline', 'pc');
end
if mode ~= 1 && mode ~= 90 && mode ~=0 && mode ~= 101 && mode ~=1000
    if issorted(ResultsMatrix(10,:)) || issorted(flip(ResultsMatrix(10,:)))
        Regression = RegressionList'\ResultsMatrix(10,:)';
        dlmwrite(strcat('LatexInputs.txt'),strcat('&',num2str(round(Regression,2))), '-append' , 'delimiter','','newline', 'pc')
    else
        dlmwrite(strcat('LatexInputs.txt'),'& -', '-append' , 'delimiter','','newline', 'pc')
    end
end
dlmwrite(strcat('LatexInputs.txt'),'\\', '-append' , 'delimiter','','newline', 'pc')

%%
dlmwrite(strcat('LatexInputs.txt'),' \textbf{2$^{nd}$ Stage Distance Flown (km)}', '-append' , 'delimiter','','newline', 'pc')
for i = 1:length(namelist)
dlmwrite(strcat('LatexInputs.txt'),strcat('& \SecondDist', namelist{i},returnparam), '-append','delimiter','','newline', 'pc');
end
if mode ~= 1 && mode ~= 90 && mode ~=0 && mode ~= 101 && mode ~=1000
    if issorted(ResultsMatrix(25,:)) || issorted(flip(ResultsMatrix(25,:)))
        Regression = RegressionList'\ResultsMatrix(25,:)';
        dlmwrite(strcat('LatexInputs.txt'),strcat('&',num2str(round(Regression,2))), '-append' , 'delimiter','','newline', 'pc')
    else
        dlmwrite(strcat('LatexInputs.txt'),'& -', '-append' , 'delimiter','','newline', 'pc')
    end
end
dlmwrite(strcat('LatexInputs.txt'),'\\', '-append' , 'delimiter','','newline', 'pc')

%%
if auxdata.returnMode == 1 % put return fuel usage into table
  
dlmwrite(strcat('LatexInputs.txt'),' \textbf{2$^{nd}$ Stage Return Fuel (kg)}', '-append' , 'delimiter','','newline', 'pc')
for i = 1:length(namelist)
dlmwrite(strcat('LatexInputs.txt'),strcat('& \returnFuel', namelist{i},returnparam), '-append','delimiter','','newline', 'pc');
end
if mode ~= 1 && mode ~= 90 && mode ~=0 && mode ~= 101 && mode ~=1000
    if issorted(ResultsMatrix(15,:)) || issorted(flip(ResultsMatrix(15,:)))
        Regression = RegressionList'\ResultsMatrix(15,:)';
        dlmwrite(strcat('LatexInputs.txt'),strcat('&',num2str(round(Regression,2))), '-append' , 'delimiter','','newline', 'pc')
    else
        dlmwrite(strcat('LatexInputs.txt'),'& -', '-append' , 'delimiter','','newline', 'pc')
    end
end
dlmwrite(strcat('LatexInputs.txt'),'\\', '-append' , 'delimiter','','newline', 'pc')  


%%
dlmwrite(strcat('LatexInputs.txt'),' \textbf{2$^{nd}$ Stage Return Distance (km)}', '-append' , 'delimiter','','newline', 'pc')
for i = 1:length(namelist)
dlmwrite(strcat('LatexInputs.txt'),strcat('& \returnDist', namelist{i},returnparam), '-append','delimiter','','newline', 'pc');
end
if mode ~= 1 && mode ~= 90 && mode ~=0 && mode ~= 101 && mode ~=1000
    if issorted(ResultsMatrix(26,:)) || issorted(flip(ResultsMatrix(26,:)))
        Regression = RegressionList'\ResultsMatrix(26,:)';
        dlmwrite(strcat('LatexInputs.txt'),strcat('&',num2str(round(Regression,2))), '-append' , 'delimiter','','newline', 'pc')
    else
        dlmwrite(strcat('LatexInputs.txt'),'& -', '-append' , 'delimiter','','newline', 'pc')
    end
end
dlmwrite(strcat('LatexInputs.txt'),'\\', '-append' , 'delimiter','','newline', 'pc')

end
dlmwrite(strcat('LatexInputs.txt'),'\hline ', '-append' , 'delimiter','','newline', 'pc')
	
%%
dlmwrite(strcat('LatexInputs.txt'),' \textbf{3$^{rd}$ Stage $\eta_{exergy}$ (\%)}', '-append' , 'delimiter','','newline', 'pc')
for i = 1:length(namelist)
dlmwrite(strcat('LatexInputs.txt'),strcat('& \textbf{\thirddExergyEff', namelist{i},returnparam,'}' ), '-append','delimiter','','newline', 'pc')
end
if mode ~= 1 && mode ~= 90 && mode ~=0 && mode ~= 101 && mode ~=1000
    if issorted(ResultsMatrix(21,:)) || issorted(flip(ResultsMatrix(21,:)))
        Regression = RegressionList'\ResultsMatrix(21,:)';
        dlmwrite(strcat('LatexInputs.txt'),strcat('& \textbf{',num2str(round(Regression,3)),'}'), '-append' , 'delimiter','','newline', 'pc')
    else
        dlmwrite(strcat('LatexInputs.txt'),'& -', '-append' , 'delimiter','','newline', 'pc')
    end
end
dlmwrite(strcat('LatexInputs.txt'),'\\', '-append' , 'delimiter','','newline', 'pc')

% %%
% dlmwrite(strcat('LatexInputs.txt'),' \textbf{3$^{rd}$ Stage Destroyed Exergy (GJ)}', '-append' , 'delimiter','','newline', 'pc')
% for i = 1:length(namelist)
% dlmwrite(strcat('LatexInputs.txt'),strcat('& \textbf{\thirdExergyDestroyed', namelist{i},returnparam,'}' ), '-append','delimiter','','newline', 'pc')
% end
% if mode ~= 1 && mode ~= 90 && mode ~=0
%     if issorted(ResultsMatrix(29,:)) || issorted(flip(ResultsMatrix(29,:)))
%         Regression = RegressionList'\ResultsMatrix(29,:)';
%         dlmwrite(strcat('LatexInputs.txt'),strcat('& \textbf{',num2str(round(Regression/10^9,3)),'}'), '-append' , 'delimiter','','newline', 'pc')
%     else
%         dlmwrite(strcat('LatexInputs.txt'),'& -', '-append' , 'delimiter','','newline', 'pc')
%     end
% end
% dlmwrite(strcat('LatexInputs.txt'),'\\', '-append' , 'delimiter','','newline', 'pc')

%%
dlmwrite(strcat('LatexInputs.txt'),' \textbf{3$^{rd}$ Stage $t$, $q >$ 5kpa (s)}', '-append' , 'delimiter','','newline', 'pc')
for i = 1:length(namelist)
dlmwrite(strcat('LatexInputs.txt'),strcat('& \thirdqOverFive', namelist{i},returnparam), '-append','delimiter','','newline', 'pc');
end
if mode ~= 1 && mode ~= 90 && mode ~=0 && mode ~= 101 && mode ~=1000
    if issorted(ResultsMatrix(11,:)) || issorted(flip(ResultsMatrix(11,:)))
        Regression = RegressionList'\ResultsMatrix(11,:)';
        dlmwrite(strcat('LatexInputs.txt'),strcat('&',num2str(round(Regression,2))), '-append' , 'delimiter','','newline', 'pc')
    else
        dlmwrite(strcat('LatexInputs.txt'),'& -', '-append' , 'delimiter','','newline', 'pc')
    end
end
dlmwrite(strcat('LatexInputs.txt'),'\\', '-append' , 'delimiter','','newline', 'pc')
%%
dlmwrite(strcat('LatexInputs.txt'),' \textbf{3$^{rd}$ Stage max $\alpha$ (deg)}', '-append' , 'delimiter','','newline', 'pc')
for i = 1:length(namelist)
dlmwrite(strcat('LatexInputs.txt'),strcat('& \thirdmaxAoA', namelist{i},returnparam), '-append','delimiter','','newline', 'pc');
end
if mode ~= 1 && mode ~= 90 && mode ~=0 && mode ~= 101 && mode ~=1000
    if issorted(ResultsMatrix(12,:)) || issorted(flip(ResultsMatrix(12,:)))
        Regression = RegressionList'\ResultsMatrix(12,:)';
        dlmwrite(strcat('LatexInputs.txt'),strcat('&',num2str(round(Regression,2))), '-append' , 'delimiter','','newline', 'pc')
    else
        dlmwrite(strcat('LatexInputs.txt'),'& -', '-append' , 'delimiter','','newline', 'pc')
    end
end
dlmwrite(strcat('LatexInputs.txt'),'\\', '-append' , 'delimiter','','newline', 'pc')

%%
dlmwrite(strcat('LatexInputs.txt'),'\textbf{3$^{rd}$ Stage Fuel Mass (kg)}', '-append' , 'delimiter','','newline', 'pc')
for i = 1:length(namelist)
dlmwrite(strcat('LatexInputs.txt'),strcat('& \thirdmFuel', namelist{i},returnparam), '-append','delimiter','','newline', 'pc');
end
if mode ~= 1 && mode ~= 90 && mode ~=0 && mode ~= 101 && mode ~=1000
    if issorted(ResultsMatrix(22,:)) || issorted(flip(ResultsMatrix(22,:)))
        Regression = RegressionList'\ResultsMatrix(22,:)';
        dlmwrite(strcat('LatexInputs.txt'),strcat('&',num2str(round(Regression,2))), '-append' , 'delimiter','','newline', 'pc')
    else
        dlmwrite(strcat('LatexInputs.txt'),'& -', '-append' , 'delimiter','','newline', 'pc')
    end
end
dlmwrite(strcat('LatexInputs.txt'),'\\', '-append' , 'delimiter','','newline', 'pc')



	dlmwrite(strcat('LatexInputs.txt'),'\hline ', '-append' , 'delimiter','','newline', 'pc')
	
dlmwrite(strcat('LatexInputs.txt'),'\end{tabular} ', '-append' , 'delimiter','','newline', 'pc')

movefile(strcat('LatexInputs.txt'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));



%% Plot Comparison Figures

if length(output)>1 && mode ~= 99
    
   
time1nom = output{ceil(length(output)/2)}.result.solution.phase(1).time.';
% 
alt1nom = output{ceil(length(output)/2)}.result.solution.phase(1).state(:,1).';
v1nom = output{ceil(length(output)/2)}.result.solution.phase(1).state(:,2).';
m1nom = output{ceil(length(output)/2)}.result.solution.phase(1).state(:,3).';
gamma1nom = output{ceil(length(output)/2)}.result.solution.phase(1).state(:,4).';
alpha1nom = output{ceil(length(output)/2)}.result.solution.phase(1).state(:,5).';
zeta1nom = output{ceil(length(output)/2)}.result.solution.phase(1).state(:,6).';
lat1nom = output{ceil(length(output)/2)}.result.solution.phase(1).state(:,8).';
lon1nom = output{ceil(length(output)/2)}.result.solution.phase(1).state(:,9).';
% 

alt21nom = output{ceil(length(output)/2)}.result.solution.phase(2).state(:,1);
alt22nom = output{ceil(length(output)/2)}.result.solution.phase(3).state(:,1);
lon21nom = output{ceil(length(output)/2)}.result.solution.phase(2).state(:,2);
lon22nom = output{ceil(length(output)/2)}.result.solution.phase(3).state(:,2);
lat21nom = output{ceil(length(output)/2)}.result.solution.phase(2).state(:,3);
lat22nom = output{ceil(length(output)/2)}.result.solution.phase(3).state(:,3);
v21nom = output{ceil(length(output)/2)}.result.solution.phase(2).state(:,4); 
v22nom = output{ceil(length(output)/2)}.result.solution.phase(3).state(:,4); 
gamma21nom = output{ceil(length(output)/2)}.result.solution.phase(2).state(:,5); 
gamma22nom = output{ceil(length(output)/2)}.result.solution.phase(3).state(:,5); 
zeta21nom = output{ceil(length(output)/2)}.result.solution.phase(2).state(:,6);
zeta22nom = output{ceil(length(output)/2)}.result.solution.phase(3).state(:,6);
alpha21nom = output{ceil(length(output)/2)}.result.solution.phase(2).state(:,7);
alpha22nom = output{ceil(length(output)/2)}.result.solution.phase(3).state(:,7);
eta21nom = output{ceil(length(output)/2)}.result.solution.phase(2).state(:,8);
eta22nom = output{ceil(length(output)/2)}.result.solution.phase(3).state(:,8);
mFuel21nom = output{ceil(length(output)/2)}.result.solution.phase(2).state(:,9); 
mFuel22nom = output{ceil(length(output)/2)}.result.solution.phase(3).state(:,9); 

throttle22nom = output{ceil(length(output)/2)}.result.solution.phase(3).state(:,10);

aoadot21nom  = output{ceil(length(output)/2)}.result.solution.phase(2).control(:,1); 
etadot21nom  = output{ceil(length(output)/2)}.result.solution.phase(2).control(:,2); 

aoadot22nom  = output{ceil(length(output)/2)}.result.solution.phase(3).control(:,1); 
etadot22nom  = output{ceil(length(output)/2)}.result.solution.phase(3).control(:,2); 

time21nom = output{ceil(length(output)/2)}.result.solution.phase(2).time;
time22nom = output{ceil(length(output)/2)}.result.solution.phase(3).time;


alt3nom  = output{ceil(length(output)/2)}.result.solution.phase(4).state(:,1);
v3nom    = output{ceil(length(output)/2)}.result.solution.phase(4).state(:,2);
gamma3nom  = output{ceil(length(output)/2)}.result.solution.phase(4).state(:,3);
m3nom    = output{ceil(length(output)/2)}.result.solution.phase(4).state(:,4);
aoa3nom    = output{ceil(length(output)/2)}.result.solution.phase(4).state(:,5);
lat3nom    = output{ceil(length(output)/2)}.result.solution.phase(4).state(:,6);
zeta3nom    = output{ceil(length(output)/2)}.result.solution.phase(4).state(:,7);
aoadot3nom       = output{ceil(length(output)/2)}.result.solution.phase(4).control(:,1);

time3nom = output{ceil(length(output)/2)}.result.solution.phase(4).time;
   
    
[rdot3,xidot3,phidot3,gammadot3,vdot3,zetadot3, mdot3, Vec_angle3, AoA_max3, T3, L3, D3, q3nom] = ThirdStageDynamics(alt3nom,gamma3nom,v3nom,m3nom,aoa3nom,time3nom,auxdata,aoadot3nom,lat3nom,zeta3nom);

[AltF_actual, v3F, altexonom, v3exo, timeexonom, mpayload, Alpha3, mexo,qexo,gammaexo,Dexo,zetaexo,latexo,incexo,Texo,CLexo,Lexo,incdiffexo,lonexo,dvtot3,m3_4] = ThirdStageSim(alt3nom(end),gamma3nom(end),v3nom(end), lat3nom(end),0, zeta3nom(end), m3nom(end), auxdata);


[altdot21,xidot21,phidot21,gammadot21,a21,zetadot21, q21nom, M21, Fd21, rho021,L21,Fueldt21,T21,Isp21,q121,flapdeflection21,heating_rate21,CG21,T121,P121,M121,P021,T021,P_1_tip21,T_1_tip21,rho_1_tip21,M_1_tip21] = SPARTANDynamics(gamma21nom, alt21nom, v21nom,auxdata,zeta21nom,lat21nom,lon21nom,alpha21nom,eta21nom,ones(length(alt21nom),1), mFuel21nom,mFuel21nom(1),mFuel21nom(end), 1, 0);
[~,~,~,~,~,~, q22nom, M22nom, Fd22, rho022,L22,Fueldt22,T22,Isp22,q122nom,flapdeflection22,heating_rate22,CG22,T122,P122,M122,P022,T022,P_1_tip22,T_1_tip22,rho_1_tip22,M_1_tip22] = SPARTANDynamics(gamma22nom, alt22nom, v22nom,auxdata,zeta22nom,lat22nom,lon22nom,alpha22nom,eta22nom,throttle22nom, mFuel22nom,0,0, 0, 0);
    
throttle22nom(M22nom<5.0) = throttle22nom(M22nom<5.0).*gaussmf(M22nom(M22nom<5.0),[.01,5]); % remove nonsense throttle points
throttle22nom(q122nom<20000) = throttle22nom(q122nom<20000).*gaussmf(q122nom(q122nom<20000),[100,20000]); % rapidly reduce throttle to 0 after passing the lower limit of 20kPa dynamic pressure. This dynamic pressure is after the conical shock.
%     





    figure(1110)
    fig = gcf;
set(fig,'Position',[200 0 850 800])
    hold on

figure(3110)
    fig = gcf;
set(fig,'Position',[200 0 850 800])
    hold on
    
    figure(2110)
    fig = gcf;
set(fig,'Position',[200 0 850 800])
    hold on
    
    
   
    
    
   figure(2210)
    fig = gcf;
set(fig,'Position',[200 0 850 800])
    hold on
%     Colours = [flip(0.7:0.3/(floor(length(output)/2)-1):1)'  zeros(floor(length(output)/2),1) zeros(floor(length(output)/2),1); 0 0 0; zeros(floor(length(output)/2),1)  zeros(floor(length(output)/2),1) (0.7:0.3/(floor(length(output)/2)-1):1)'];
%     
    Colours = [0.8500, 0.3250, 0.0980; 1 0 0; 0 0 0; 0 0 1; 0.3010, 0.7450, 0.9330];
    
    LineStyleList = {'--','-','-','-','--'};
    
    for j = 1:length(output)
        
time1 = output{j}.result.solution.phase(1).time.';
% 
alt1 = output{j}.result.solution.phase(1).state(:,1).';
v1 = output{j}.result.solution.phase(1).state(:,2).';
m1 = output{j}.result.solution.phase(1).state(:,3).';
gamma1 = output{j}.result.solution.phase(1).state(:,4).';
alpha1 = output{j}.result.solution.phase(1).state(:,5).';
zeta1 = output{j}.result.solution.phase(1).state(:,6).';
lat1 = output{j}.result.solution.phase(1).state(:,8).';
lon1 = output{j}.result.solution.phase(1).state(:,9).';
% 

alt21 = output{j}.result.solution.phase(2).state(:,1);
alt22 = output{j}.result.solution.phase(3).state(:,1);
lon21 = output{j}.result.solution.phase(2).state(:,2);
lon22 = output{j}.result.solution.phase(3).state(:,2);
lat21 = output{j}.result.solution.phase(2).state(:,3);
lat22 = output{j}.result.solution.phase(3).state(:,3);
v21 = output{j}.result.solution.phase(2).state(:,4); 
v22 = output{j}.result.solution.phase(3).state(:,4); 
gamma21 = output{j}.result.solution.phase(2).state(:,5); 
gamma22 = output{j}.result.solution.phase(3).state(:,5); 
zeta21 = output{j}.result.solution.phase(2).state(:,6);
zeta22 = output{j}.result.solution.phase(3).state(:,6);
alpha21 = output{j}.result.solution.phase(2).state(:,7);
alpha22 = output{j}.result.solution.phase(3).state(:,7);
eta21 = output{j}.result.solution.phase(2).state(:,8);
eta22 = output{j}.result.solution.phase(3).state(:,8);
mFuel21 = output{j}.result.solution.phase(2).state(:,9); 
mFuel22 = output{j}.result.solution.phase(3).state(:,9); 

throttle22 = output{j}.result.solution.phase(3).state(:,10);



aoadot21  = output{j}.result.solution.phase(2).control(:,1); 
etadot21  = output{j}.result.solution.phase(2).control(:,2); 

aoadot22  = output{j}.result.solution.phase(3).control(:,1); 
etadot22  = output{j}.result.solution.phase(3).control(:,2); 

time21 = output{j}.result.solution.phase(2).time;
time22 = output{j}.result.solution.phase(3).time;


alt3  = output{j}.result.solution.phase(4).state(:,1);
v3    = output{j}.result.solution.phase(4).state(:,2);
gamma3  = output{j}.result.solution.phase(4).state(:,3);
m3    = output{j}.result.solution.phase(4).state(:,4);
aoa3    = output{j}.result.solution.phase(4).state(:,5);
lat3    = output{j}.result.solution.phase(4).state(:,6);
zeta3    = output{j}.result.solution.phase(4).state(:,7);
aoadot3       = output{j}.result.solution.phase(4).control(:,1);

time3 = output{j}.result.solution.phase(4).time;

        
[rdot3,xidot3,phidot3,gammadot3,vdot3,zetadot3, mdot3, Vec_angle3, AoA_max3, T3, L3, D3, q3] = ThirdStageDynamics(alt3,gamma3,v3,m3,aoa3,time3,auxdata,aoadot3,lat3,zeta3);

lon3 = [];
lon3(1) = lon21(end);
for i = 2:length(time3)
    lon3(i) = lon3(i-1) + xidot3(i-1)*(time3(i)-time3(i-1));
end


[altdot21,xidot21,phidot21,gammadot21,a21,zetadot21, q21, M21, Fd21, rho021,L21,Fueldt21,T21,Isp21,q121,flapdeflection21,heating_rate21,CG21] = SPARTANDynamics(gamma21, alt21, v21,auxdata,zeta21,lat21,lon21,alpha21,eta21,ones(length(alt21),1), mFuel21,mFuel21(1),mFuel21(end), 1, 0);
[~,~,~,~,~,~, q22, M22, Fd22, rho022,L22,Fueldt22,T22,Isp22,q122,flapdeflection22,heating_rate22] = SPARTANDynamics(gamma22, alt22, v22,auxdata,zeta22,lat22,lon22,alpha22,eta22,throttle22, mFuel22,0,0, 0, 0);
% 
[AltF_actual, v3F, altexo, v3exo, timeexo, mpayload, Alpha3, mexo,qexo,gammaexo,Dexo,zetaexo,latexo,incexo,Texo,CLexo,Lexo,incdiffexo,lonexo] = ThirdStageSim(alt3(end),gamma3(end),v3(end), lat3(end),lon3(end), zeta3(end), m3(end), auxdata);

throttle22(M22<5.0) = throttle22(M22<5.0).*gaussmf(M22(M22<5.0),[.01,5]); % remove nonsense throttle points
throttle22(q122<20000) = throttle22(q122<20000).*gaussmf(q122(q122<20000),[100,20000]); % rapidly reduce throttle to 0 after passing the lower limit of 20kPa dynamic pressure. This dynamic pressure is after the conical shock.
%     




        figure(3110)
    subplot(4,1,1)
    hold on
    title('Altitude (km','FontSize',9);
    plot([time3-time3(1); timeexo.'+time3(end)-time3(1)], [alt3; altexo.']/1000','Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])
%     xlim([time3(1) timeexo(end)+time3(end)])

% subplot(4,2,2)
%     hold on
%     title('Altitude Difference (km)','FontSize',9);
% 
%     plot(0:1/(length(time3)+length(timeexo)-1):1, ([alt3; altexo.']'-interp1(0:1/(length(time3nom)+length(timeexonom)-1):1,[alt3nom; altexonom.'],0:1/(length(time3)+length(timeexo)-1):1))/1000','Color',Colours(j,:),'LineStyle',LineStyleList{j})
% set(gca,'xticklabels',[])
% %     xlim([time3(1) timeexo(end)+time3(end)])


    subplot(4,1,3)
    hold on
    title('Angle of Attack (deg)','FontSize',9);
    plot([time3-time3(1); timeexo.'+time3(end)-time3(1)],[rad2deg(aoa3);0*ones(length(timeexo),1)],'Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])

% subplot(4,2,4)
%     hold on
%     title('Angle of Attack Difference (deg)','FontSize',9);
%     plot(0:1/(length(time3)-1):1, rad2deg(aoa3'-interp1(0:1/(length(time3nom)-1):1,aoa3nom,0:1/(length(time3)-1):1))','Color',Colours(j,:),'LineStyle',LineStyleList{j})
% set(gca,'xticklabels',[])

%     xlim([time3(1) timeexo(end)+time3(end)])
    subplot(4,1,2)
    hold on
    title('Velocity (m/s)','FontSize',9);
    plot([time3-time3(1); timeexo.'+time3(end)-time3(1)],[v3;v3exo.'],'Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])

% subplot(4,2,6)
%     hold on
%     title('Velocity Difference (m/s)','FontSize',9);
%     plot(0:1/(length(time3)-1):1, (v3'-interp1(0:1/(length(time3nom)-1):1,v3nom,0:1/(length(time3)-1):1))','Color',Colours(j,:),'LineStyle',LineStyleList{j})
% set(gca,'xticklabels',[])


%     xlabel('Time (s)','FontSize',9);
%     xlim([time3(1) timeexo(end)+time3(end)])
    subplot(4,1,4)
    hold on
    title('Trajectory Angle (deg)','FontSize',9);
    plot([time3-time3(1); timeexo.'+time3(end)-time3(1)], [rad2deg(gamma3);rad2deg(gammaexo).'],'Color',Colours(j,:),'LineStyle',LineStyleList{j})

    xlabel('Time (s)','FontSize',9);
    
    
%     subplot(4,2,8)
%     hold on
%     title('Flight Path Angle Difference (Deg)','FontSize',9);
%     plot(0:1/(length(time3)-1):1, rad2deg(gamma3'-interp1(0:1/(length(time3nom)-1):1,gamma3nom,0:1/(length(time3)-1):1))','Color',Colours(j,:),'LineStyle',LineStyleList{j})
% % set(gca,'xticklabels',[])
% xlabel('Normalised Time','FontSize',9);
% %     xlim([time3(1) timeexo(end)+time3(end)])



        figure(2110)
    subplot(4,1,1)
    hold on
    title('Altitude (km','FontSize',9);
    plot([time21-time21(1)], [alt21]/1000','Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])
%     xlim([time3(1) timeexo(end)+time3(end)])

% subplot(3,2,2)
%     hold on
%     title('Altitude Difference (km)','FontSize',9);

 
    plot(0:1/(length(time21)-1):1, (alt21'-interp1(0:1/(length(time21nom)-1):1,alt21nom,0:1/(length(time21)-1):1))/1000','Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])
%     xlim([time3(1) timeexo(end)+time3(end)])


    subplot(4,1,2)
    hold on
    title('Velocity (m/s)','FontSize',9);
    plot([time21-time21(1)],[v21],'Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])

    subplot(4,1,3)
    hold on
    title('Dynamic Pressure (kPa)','FontSize',9);
    plot([time21-time21(1)],[q21]/1000,'Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])

  subplot(4,1,4)
    hold on
    title('Angle of Attack (deg)','FontSize',9);
    plot([time21-time21(1)],rad2deg(alpha21),'Color',Colours(j,:),'LineStyle',LineStyleList{j})
xlabel('Time (s)');
    
    
%     xlim([time21(1) timeexo(end)+time21(end)])

% subplot(3,2,4)
%     hold on
%     title('Dynamic Pressure Difference (kPa)','FontSize',9);
%     plot(0:1/(length(time21)-1):1, (q21'-interp1(0:1/(length(time21nom)-1):1,q21nom,0:1/(length(time21)-1):1))/1000','Color',Colours(j,:),'LineStyle',LineStyleList{j})
% set(gca,'xticklabels',[])





% set(gca,'xticklabels',[])


% subplot(3,2,6)
%     hold on
%     title('Velocity Difference (m/s)','FontSize',9);
%     plot(0:1/(length(time21)-1):1, (v21'-interp1(0:1/(length(time21nom)-1):1,v21nom,0:1/(length(time21)-1):1))','Color',Colours(j,:),'LineStyle',LineStyleList{j})
% set(gca,'xticklabels',[])


        figure(2210)
    subplot(4,1,1)
    hold on
    title('Altitude (km','FontSize',9);
    plot([time22-time22(1)], [alt22]/1000','Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])
%     xlim([time3(1) timeexo(end)+time3(end)])

% subplot(4,2,2)
%     hold on
%     title('Altitude Difference (km)','FontSize',9);

 
%     plot(time22-time22(1), (alt22'-interp1(time22nom,alt22nom,time22)')/1000,'Color',Colours(j,:),'LineStyle',LineStyleList{j})
% set(gca,'xticklabels',[])
%     xlim([time3(1) timeexo(end)+time3(end)])


    subplot(4,1,2)
    hold on
    title('Dynamic Pressure (kPa','FontSize',9);
    plot([time22-time22(1)],[q22]/1000,'Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])
%     xlim([time22(1) timeexo(end)+time22(end)])

% subplot(4,2,4)
%     hold on
%     title('Dynamic Pressure Difference (kPa)','FontSize',9);
%     plot(time22-time22(1), (q22'-interp1(time22nom,q22nom,time22)')/1000,'Color',Colours(j,:),'LineStyle',LineStyleList{j})
% set(gca,'xticklabels',[])


%     xlim([time22(1) timeexo(end)+time22(end)])
    subplot(4,1,3)
    hold on
    title('Velocity (m/s)','FontSize',9);
    plot([time22-time22(1)],[v22],'Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])

% xlabel('Time (s)','FontSize',9);

% subplot(4,2,6)
%     hold on
%     title('Velocity Difference (m/s)','FontSize',9);
%     plot(time22-time22(1), (v22'-interp1(time22nom,v22nom,time22)'),'Color',Colours(j,:),'LineStyle',LineStyleList{j})
% set(gca,'xticklabels',[])

% xlabel('Normalised Time','FontSize',9);


subplot(4,1,4)
    hold on
    title('Throttle','FontSize',9);
    plot([time22-time22(1)], [throttle22],'Color',Colours(j,:),'LineStyle',LineStyleList{j})

    xlabel('Time (s)','FontSize',9);

    end
    figure(2110)
    legend(namelist)
    figure(2210)
    legend(namelist)
    figure(3110)
    legend(namelist)
    
saveas(figure(2110),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,'SecondStageComparison.fig']);
print(figure(2110),'SecondStageComparison','-dpng');
movefile('SecondStageComparison.png',sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));

saveas(figure(2210),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,'ReturnComparison.fig']);
print(figure(2210),'ReturnComparison','-dpng');
movefile('ReturnComparison.png',sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));


saveas(figure(3110),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))),filesep,'ThirdStageComparison.fig']);
print(figure(3110),'ThirdStageComparison','-dpng');
movefile('ThirdStageComparison.png',sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));


end
close all




if mode == 99
    for i = 1:length(output)
    payload_list(i) = -output{i}.result.objective;
    
    end
    payload_list = payload_list';
    
    
    figure(9999)
    
    interactionplot(payload_list,auxdata.p,'varnames',{'q_max','Isp','Cd'})
    
end



%% Plot Mesh

if mode == 1
    
    figure(9001);
    
    
for j = 1:4
    
    subplot(2,2,j)
    
    
    solution = output{1}.result.solution;
time = solution.phase(j).time;
state = solution.phase(j).state;
control = solution.phase(j).control;
for i=1:length(output{1}.meshhistory);
  mesh(i).meshPoints = [0 cumsum(output{1}.meshhistory(i).result.setup.mesh.phase(j).fraction)];
  mesh(i).time =  output{1}.meshhistory(i).result.solution.phase(j).time;
  mesh(i).iteration = i*ones(size(mesh(i).meshPoints));
  mesh(i).iterationTime = i*ones(size(mesh(i).time));
end;

for i=1:length(mesh);
  pp = plot(mesh(i).meshPoints,mesh(i).iteration,'bo');
  set(pp,'LineWidth',1.25);
  hold on;
end;
xl = xlabel('Mesh Point Location');
yl = ylabel('Mesh Iteration');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'YTick',0:length(mesh),'FontSize',16,'FontName','Times');
if auxdata.returnMode == 1
 axis([0 1 1 5]);   
else
axis([0 1 1 4]);
end
grid on;
print -dpng MeshHistory.png
movefile('MeshHistory.png',sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode),num2str(returnMode))));

end
end

close all
end




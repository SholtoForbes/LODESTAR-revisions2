function [Kn,Fq_axial_rarified] = ThirdStageHAViscous(T0,P0,M0,alpha)
% This function calculates the viscous effects on the third stage rocket at
% high altitude, where the flow becomes rarified. From Multi-Disciplinary Modelling of
% Future Space-Access Vehicles by Romain Wuilbercq

% It is assumed that the whole nozzle is shadowed

cone_angle = deg2rad(10.3885);
d = 1.1;
nose_rad = 0.0712;
L= 8.66;

K_B = 1.380648813e-23; %Boltzmann constant
d_air = 4e-10; % Diameter of dry air

Kn = K_B.*T0./(sqrt(2).*pi.*d_air.^2.*P0.*L); % Calculate Knudsen no.

gamma = 1.4;
sigma_S = sqrt(gamma/2).*M0;


% Icosagonal Nose Tip
panel_no = 1;

R_icon = 0.001; % radius of icosagonal nose tip
inc(panel_no,:) = deg2rad(90) + alpha;
A(panel_no) = 5*R_icon^2/2*(sqrt(5)-1);

%Nose
no_partitions_nose = 10; % Number of lengthwise partitions
end_angle_hemisphere = deg2rad(90) - cone_angle; %angle that the spherical tip ends at

for i = 1:20 %20 because of icosagonal tip, starts at panel to one side of top
    for j = 1:no_partitions_nose
        inc(panel_no,:) = end_angle_hemisphere/(no_partitions_nose+1)*j + alpha*cos(2*pi/20*i);  % calculate inclinations for each panel, not including nose, which counts at 0 degrees
        
        l_forward = 2.*pi.*0.0711667.*sin(end_angle_hemisphere/(no_partitions_nose+1)*j)/20 ; % length of foward panel side, appriximated as circular circumference. 0.0711667 is radius of full sphere (which the nose is not)
        l_rear = 2.*pi.*0.0711667.*sin(end_angle_hemisphere/(no_partitions_nose+1)*(j+1))/20 ; % length of rear panel side
        l_side = (0.0711667.*(pi/2 - cone_angle) - R_icon) ./ no_partitions_nose; % length of side approximated as part of a cirucumference. quarter sphere minus the cone angle and the tip
        A(panel_no) = (l_forward + l_rear)/2.*l_side;
        
        panel_no = panel_no + 1;
    end
end

% Cone

for i = 1:20
    inc(panel_no,:) = cone_angle + alpha*cos(2*pi/20*i);
    A(panel_no) = (2.*pi.*0.0711667.*sin(end_angle_hemisphere) + 2.*pi.*0.55)/2.*2.6618/20; % average of circumference at start and end of cone. length of 2.6618 measured from creo
    panel_no = panel_no + 1;
end

%Cylindrical Body
for i = 1:20
    inc(panel_no,:) = alpha*cos(2*pi/20*i);
    A(panel_no) = 2.*pi.*0.55.*4.5/20;
    panel_no = panel_no + 1;
end

% Calculate skin friction coefficients from inclinations
for i = 1:length(inc(:,1))
    
    Cf_rar(i,:) = cos(inc(i,:))./(2.*sqrt(pi).*sigma_S).*(exp(-(sigma_S.*sin(inc(i,:))).^2)+sqrt(pi).*(sigma_S.*sin(inc(i,:))).*(1+erf(sigma_S.*sin(inc(i,:)))));
    % Diffuse reflection is assumed
    Shear = Cf_rar;
    
end


for j = 1:length(Shear(1,:))
    Fq_axial_rarified(j) = Shear(:,j)'*A';
end

end


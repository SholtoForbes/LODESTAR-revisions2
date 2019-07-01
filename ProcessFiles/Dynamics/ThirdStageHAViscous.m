function [Kn] = ThirdStageHAViscous(T0,P0,L)
% This function calculates the viscous effects on the third stage rocket at
% high altitude, where the flow becomes rarified. From Multi-Disciplinary Modelling of
% Future Space-Access Vehicles by Romain Wuilbercq

K_B = 1.380648813e-23; %Boltzmann constant
d_air = 4e-10; % Diameter of dry air

Kn = K_B.*T0./(sqrt(2).*pi.*d_air.^2.*P0.*L); % Calculate Knudsen no.

end


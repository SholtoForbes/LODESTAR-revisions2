function [rdot,xidot,phidot,gammadot,vdot,zetadot] = RotCoordsRocket(alt,xi,phi,gamma,v,zeta,L,D,T,m,alpha,Vec_angle)
% % Determination of motion in rotating coordinates
% %xi  Longitude (rad)
% %phi  Latitude (rad)
% %gamma  Flight Path Angle (rad)
% %zeta  Heading Angle (rad)
% mu_E = 3.986e14; % m^3/s^2 Earth Gravitational Parameter
% omega_E = 7.292115e-5; % s^-1 Earth Rotation Rate
% 
% rdot = v.*sin(gamma);
% 
% xidot = v.*cos(gamma).*cos(zeta)./(r.*cos(phi));
% 
% phidot = v.*cos(gamma).*sin(zeta)./r;
% 
% gammadot = T.*sin(alpha+Vec_angle)./(m.*v) + (v./r - mu_E./(r.^2.*v)).*cos(gamma) + L./(m.*v) + cos(phi).*(2.*omega_E.*cos(zeta) + omega_E.^2.*r./v.*(cos(phi).*cos(gamma)+sin(phi).*sin(gamma).*sin(zeta)));
% 
% vdot = T.*cos(alpha+Vec_angle)./(m) - mu_E.*sin(gamma)./r.^2 -D./m + omega_E.^2.*r.*cos(phi).*(cos(phi).*cos(gamma)+sin(phi).*sin(gamma).*sin(zeta));
% 
% zetadot = -v./r.*tan(phi).*cos(gamma).*cos(zeta) + 2.*omega_E.*cos(phi).*tan(gamma).*sin(zeta) - omega_E.^2.*r./(v.*cos(gamma)).*sin(phi).*cos(phi).*cos(zeta)-2.*omega_E.*sin(phi);

eta = 0;


omega_E = 7.292115e-5; % s^-1 Earth Rotation Rate


Re = geocradius(rad2deg(phi)); %Calculate 


r = alt+Re;


phi_geod = geoc2geod(rad2deg(phi), r); % calculate geodetic latitude

[gn, gt] = gravitywgs84( alt, phi_geod, rad2deg(xi), 'Exact', 'Warning'); % calculate normal and tangential components of gravity
% gn = gravitywgs84( alt, phi_geod, 'TaylorSeries', 'None'); % calculate normal and tangential components of gravity

% gt = 0;

% I may need to change this to calculate in a different way?
% potentially gravitysphericalharmonic with ecef2aer or ecef2enuv

% gn = 9.7;
% gt = 0;

rdot = v.*sin(gamma);


xidot = v.*cos(gamma).*cos(zeta)./(r.*cos(phi));

phidot = v.*cos(gamma).*sin(zeta)./r;

total_lift = T.*sin(alpha+Vec_angle) + L;
% total_lift = L;

gammadot = total_lift./(m.*v).*cos(eta) + (v./r - gn./v).*cos(gamma) + cos(phi).*(2.*omega_E.*cos(zeta) + omega_E.^2.*r./v.*(cos(phi).*cos(gamma)+sin(phi).*sin(gamma).*sin(zeta))) - gt./v.*sin(gamma).*cos(zeta);


vdot = T.*cos(alpha+Vec_angle)./(m) - gn.*sin(gamma) -D./m + omega_E.^2.*r.*cos(phi).*(sin(gamma).*cos(phi) - cos(gamma).*sin(zeta).*sin(phi)) + gt.*cos(gamma).*cos(zeta);


zetadot = total_lift./(m.*v).*sin(eta)./cos(gamma) -v./r.*tan(phi).*cos(gamma).*cos(zeta) + 2.*omega_E.*cos(phi).*tan(gamma).*sin(zeta) - omega_E.^2.*r./(v.*cos(gamma)).*sin(gamma).*cos(phi).*cos(zeta)-2.*omega_E.*sin(phi) - gt.*sin(zeta);

end


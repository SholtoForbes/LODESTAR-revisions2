function [rdot,xidot,phidot,gammadot,vdot,zetadot,total_lift] = RotCoords(alt,xi,phi,gamma,v,zeta,L,D,T,m,alpha,eta,delta,auxdata)
% Determination of motion in rotating coordinates with WGS84 correction,
% maddock 2017
% r radius from centre of Earth (m)
%xi  Longitude (rad)
%phi  Latitude (rad)
%gamma  Flight Path Angle (rad)
%zeta  Heading Angle (rad)
% alpha Angle of Attack (deg)
% eta ROll Angle (rad)
% L Lift (N)
% D Drag (N)
% m Mass (kg)


% alpha = deg2rad(alpha);

% mu_E = 3.986e14; % m^3/s^2 Earth Gravitational Parameter
omega_E = 7.292115e-5; % s^-1 Earth Rotation Rate


Re = geocradius(rad2deg(phi)); %Calculate 

r = alt+Re;

% phi_geod = geoc2geod(rad2deg(phi), r); % 
% 
% [gn, gt] = gravitywgs84( alt, phi_geod, rad2deg(xi), 'Exact', 'Warning'); % calculate normal and tangential components of gravity

gn = auxdata.interp.gn_interp(rad2deg(phi), alt);
gt = auxdata.interp.gt_interp(rad2deg(phi), alt);

rdot = v.*sin(gamma);


xidot = v.*cos(gamma).*cos(zeta)./(r.*cos(phi));

phidot = v.*cos(gamma).*sin(zeta)./r;

total_lift = T.*sin(alpha+delta) + L;
% total_lift = L;

gammadot = total_lift./(m.*v).*cos(eta) + (v./r - gn./v).*cos(gamma) + cos(phi).*(2.*omega_E.*cos(zeta) + omega_E.^2.*r./v.*(cos(phi).*cos(gamma)+sin(phi).*sin(gamma).*sin(zeta))) - gt./v.*sin(gamma).*cos(zeta);


vdot = T.*cos(alpha+delta)./(m) - gn.*sin(gamma) -D./m + omega_E.^2.*r.*cos(phi).*(sin(gamma).*cos(phi) - cos(gamma).*sin(zeta).*sin(phi)) + gt.*cos(gamma).*cos(zeta);


zetadot = total_lift./(m.*v).*sin(eta)./cos(gamma) -v./r.*tan(phi).*cos(gamma).*cos(zeta) + 2.*omega_E.*cos(phi).*tan(gamma).*sin(zeta) - omega_E.^2.*r./(v.*cos(gamma)).*sin(gamma).*cos(phi).*cos(zeta)-2.*omega_E.*sin(phi) - gt.*sin(zeta);

end


function [a, e, i, OM, om, th]=car2par(rr, vv, mu)

% Trasformation from cartesian coordinates to Keplerian parameters
%
% [a, e, i, OM, om, th]=car2par(rr, vv, mu)
%
% --------------------------------------------------------------------
% Input arguments:
% rr        [3x1]   Position vector             [km]
% vv        [3x1]   Velocity vector             [km/s]
% mu        [1x1]   Gravitational parameter     [km^3/s^s]
%
% -------------------------------------------------------------------
% Output arguments:
% a         [1x1]   Semi-major axis             [km]
% e         [1x1]   Eccentricity                [-]
% i         [1x1]   Inclination                 [rad]
% OM        [1x1]   RAAN                        [rad]
% om        [1x1]   Pericenter anomaly          [rad]
% th        [1x1]   True anomaly                [rad]
%
% ----------------------------------------------------------------------

% 1. Moduli posizione e velocità

r = norm(rr);
v = norm(vv);

% 2. Semiasse maggiore

a = 1/(2/r - ((v^2)/mu));

% 3. Vettore momento angolare specifico

hh = cross(rr, vv);
h  = norm(hh);

% 4. Vettore eccentricità ed eccentricità

ee = cross(vv, hh)/mu - rr/r;
e  = norm(ee);

% 5. Inclinazione

i = acos(hh(3)/h);

% 6. Linea dei nodi

N = cross([0; 0; 1], hh)/norm(cross([0; 0; 1], hh));

% 7. Ascensione retta del nodo ascendente (RAAN)

if N(2)<0
    OM = 2*pi-acos(N(1));
else
    OM = acos(N(1));
end

%OM = rad2deg(OM);

% 8. Anomalia del pericentro

if ee(3)<0
    om = 2*pi - acos(dot(N,ee)/e);
else
    om = acos(dot(N,ee)/e);
end

%om = rad2deg(om);

% 9. Anomalia vera

vr = dot(vv, rr)/r;

if vr<0
    th = 2*pi - acos(dot(rr, ee)/(r*e));
else
    th = acos(dot(rr, ee)/(r*e));
end

%th = rad2deg(th);








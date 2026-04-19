function [rr, vv]=par2car(a, e, i, OM, om, th, mu)

% Trasformation from Keplerian parameters to cartesian coordinates
%
% [rr, vv] = par2car(a, e, i, OM, om, th, mu)
%
% --------------------------------------------------------------------
% Input arguments:
% a         [1x1]   Semi-major axis             [km]
% e         [1x1]   Eccentricity                [-]
% i         [1x1]   Inclination                 [rad]
% OM        [1x1]   RAAN                        [rad]
% om        [1x1]   Pericenter anomaly          [rad]
% th        [1x1]   True anomaly                [rad]
% mu        [1x1]   Gravitational parameter     [km^3/s^s]
%
% ---------------------------------------------------------------------
% Output arguments:
% rr        [3x1]   Position vector             [km]
% vv        [3x1]   Velocity vector             [km/s]
%
% ----------------------------------------------------------------------

% 1. Semilato retto

p = a*(1-e^2);

% 2. Modulo della posizione

r = p/(1 + e*cos(th));

% 3. Vettore [rt, vt]

rt = r.*[cos(th); sin(th); 0];
vt = sqrt(mu/p).*[-sin(th); e+cos(th); 0];

% 4. Vettore [rr; vv]

% Matrici di rotazione

R_OM = [ cos(OM)  sin(OM)  0;
        -sin(OM)  cos(OM)  0;
             0       0     1];

R_i = [    1     0        0   ;
           0   cos(i)  sin(i) ;
           0  -sin(i)  cos(i) ];

R_om = [ cos(om)  sin(om)  0;
        -sin(om)  cos(om)  0;
             0       0     1];

T = R_om * R_i * R_OM;

rr = T'*rt;
vv = T'*vt;











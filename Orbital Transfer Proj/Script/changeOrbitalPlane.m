function [DeltaV, omf, theta] = changeOrbitalPlane(a, e, i_i, OMi, omi, i_f, OMf, mu)

% Change of Plane maneuvre
%
% [DeltaV, omf, theta] = changeOrbitalPlane(a, e, i_i, OMi, omi, i_f, OMf, mu)
%
% --------------------------------------------------------------------
% Input arguments:
% a         [1x1]   semi-major axis                    [km]
% e         [1x1]   eccentricity                       [-]
% i_i       [1x1]   initial inclination                [km]
% OMi       [1x1]   initial RAAN                       [-]
% omi       [1x1]   initial pericenter anomaly         [km]  
% i_f       [1x1]   final inclination                  [km^3/s^s]
% OMf       [1x1]   final RAAN                         [km]  
% mu        [1x1]   Gravitational parameter            [km^3/s^s]
%
% ---------------------------------------------------------------------
% Output arguments:
% DeltaV        [1x1]   maneuver impulse                    [km/s]
% omf           [1x1]   final pericenter anomaly            [rad]
% theta         [1x1]   true anomaly at maneuver            [rad]
%
% ----------------------------------------------------------------------

dOM   = OMf - OMi;
d_i   = i_f - i_i;
alpha = acos( cos(i_i)*cos(i_f) + sin(i_i)*sin(i_f)*cos(abs(dOM)) );

if dOM > 0 && d_i > 0

    c_ui = ( -cos(i_f) + cos(alpha)*cos(i_i) ) / (sin(alpha)*sin(i_i));
    c_uf = (  cos(i_i) - cos(alpha)*cos(i_f) ) / (sin(alpha)*sin(i_f));

    s_ui = ( sin(dOM) * sin(i_f)) / sin(alpha);
    s_uf = ( sin(dOM) * sin(i_i)) / sin(alpha);

    u_i = atan2 ( s_ui, c_ui);
    u_f = atan2 ( s_uf, c_uf);

    theta = u_i - omi;
    omf   = u_f - theta; 

elseif dOM > 0 && d_i < 0

    c_ui = (   cos(i_f) - cos(alpha)*cos(i_i) ) / (sin(alpha)*sin(i_i));
    c_uf = (  -cos(i_i) + cos(alpha)*cos(i_f) ) / (sin(alpha)*sin(i_f));

    s_ui = ( sin(dOM) * sin(i_f)) / sin(alpha);
    s_uf = ( sin(dOM) * sin(i_i)) / sin(alpha);

    u_i = atan2 ( s_ui, c_ui);
    u_f = atan2 ( s_uf, c_uf);

    theta = 2*pi - u_i - omi;
    omf   = 2*pi - u_f - theta;

elseif dOM < 0 && d_i > 0

    dOM = abs(dOM);

    c_ui = (  -cos(i_f) + cos(alpha)*cos(i_i) ) / (sin(alpha)*sin(i_i));
    c_uf = (   cos(i_i) - cos(alpha)*cos(i_f) ) / (sin(alpha)*sin(i_f));

    s_ui = ( sin(abs(dOM)) * sin(i_f)) / sin(alpha);
    s_uf = ( sin(abs(dOM)) * sin(i_i)) / sin(alpha);

    u_i = atan2 ( s_ui, c_ui);
    u_f = atan2 ( s_uf, c_uf);

    theta = 2*pi - u_i - omi;
    omf   = 2*pi - u_f - theta;

elseif dOM < 0 && d_i < 0

    dOM = abs(dOM);

    c_ui = (   cos(i_f) - cos(alpha)*cos(i_i) ) / (sin(alpha)*sin(i_i));
    c_uf = (  -cos(i_i) + cos(alpha)*cos(i_f) ) / (sin(alpha)*sin(i_f));

    s_ui = ( sin(abs(dOM)) * sin(i_f)) / sin(alpha);
    s_uf = ( sin(abs(dOM)) * sin(i_i)) / sin(alpha);

    u_i = atan2 ( s_ui, c_ui);
    u_f = atan2 ( s_uf, c_uf);

    theta = u_i - omi;
    omf   = u_f - theta;

end

if cos(theta) > 0 

    theta = theta + pi;
end


V_th = ( mu / (a * (1-e^2)))^(1/2) * (1 + e*cos(theta));

DeltaV = 2 * V_th * sin(alpha/2);























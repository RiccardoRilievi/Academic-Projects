function [DeltaV, thi, thf] = changePericenterArg (a, e, omi, omf, mu)

% Change of Pericenter Argument maneuver
%
% [DeltaV, thi, thf] = changePericenterArg (a, e, omi, omf, mu)
%
% --------------------------------------------------------------------
% Input arguments:
% a         [1x1]   Semi-major axis              [km]
% e         [1x1]   Eccentricity                 [-]
% i         [1x1]   Inclination                  [rad]
% omi       [1x1]   initial pericenter anomaly   [rad]
% omf       [1x1]   final ericenter anomaly      [rad]
% th        [1x1]   True anomaly                 [rad]
% mu        [1x1]   Gravitational parameter      [km^3/s^s]
%
% ---------------------------------------------------------------------
% Output arguments:
% DeltaV    [1x1]   maneuver impulse             [km/s]
% thi       [2x1]   initial true anomaly         [rad]
% thf       [2x1]   final true anomaly           [rad]
%
% ----------------------------------------------------------------------

% Anomalia vera della manovra

d_om = omf-omi;

thi = [d_om/2; pi+d_om/2];
thf = [2*pi-d_om/2; pi-d_om/2];

for j = 1:length(thi)

    if thi(j) <0;

        thi(j) = thi(j) + 2*pi;
        
    end
    
    if thi(j) > 2*pi
        
        thi (j) = thi(j) - 2*pi;
    end

end

for j = 1:length(thf)

    if thf(j) <0;

        thf(j) = thf(j) + 2*pi;
        
    end
    
    if thf(j) > 2*pi
        
        thf(j) = thf(j) - 2*pi;
    end

end


% Costo Manovra

DeltaV = 2 * sqrt( mu / (a*(1-e^2))) * e * sin(d_om/2);

















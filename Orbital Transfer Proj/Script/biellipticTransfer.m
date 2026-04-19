function [DeltaV1, DeltaV2, DeltaV3, Deltat1, Deltat2] = biellipticTransfer(ai, ei, af, ef, ra_t, mu)

% Bitangent transfer for elliptic orbits
%
% [DeltaV1, DeltaV2, DeltaV3, Deltat1, Deltat2] = bitangetTransfer(ai, ei, af, ef, ra_t, mu)
%
% --------------------------------------------------------------------
% Input arguments:
% ai         [1x1]   initial semi-major axis            [km]
% ei         [1x1]   Initial eccentricity               [-]
% af         [1x1]   Final semi-major axis              [km]
% ef         [1x1]   Final eccentricity                 [-]
% ra_t       [1x1]   transfer orbits apocenter distance [km]  
% mu         [1x1]   Gravitational parameter            [km^3/s^s]
%
% ---------------------------------------------------------------------
% Output arguments:
% DeltaV1       [1x1]   1st maneuver impulse            [km/s]
% DeltaV2       [1x1]   2nd maneuver impulse            [km/s]
% DeltaV3       [1x1]   2rd maneuver impulse            [km/s]
% Deltat1       [1x1]   maneuvre time 1                 [s]
% Deltat2       [1x1]   maneuvre time 2                 [s]
%
% ----------------------------------------------------------------------

rp_i = ai*(1-ei);
ra_T1 = ra_t;

rp_T1 = rp_i;
ra_T2 = ra_T1;
rp_T2 = af*(1-ef);

a_T1 = (rp_T1+ra_T1)/2;
a_T2 = (rp_T2+ra_T2)/2;

e_T1 = (ra_T1-rp_T1)/(rp_T1+ra_T1);
e_T2 = (ra_T2-rp_T2)/(rp_T2+ra_T2);

DeltaV1 = mu^(1/2)*((2/rp_T1-1/a_T1)^(1/2)-(2/rp_i-1/ai)^(1/2));
DeltaV2 = mu^(1/2)*((2/ra_T2-1/a_T2)^(1/2)-(2/ra_T1-1/a_T1)^(1/2));
DeltaV3 = mu^(1/2)*((2/rp_T2-1/af)^(1/2)-(2/rp_T2-1/a_T2)^(1/2));

Deltat1  = pi*((a_T1^3)/mu)^(1/2);
Deltat2  = pi*((a_T2^3)/mu)^(1/2);




function [DeltaV1, DeltaV2, Deltat, a_t, e_t] = bitangetTransfer(a_i, e_i, a_f, e_f, type, mu)

% Bitangent transfer for elliptic orbits
%
% [DeltaV1, DeltaV2, Deltat] =  bitangetTransfer(a_i, e_i, a_f, e_f, type, mu)
%
% --------------------------------------------------------------------
% Input arguments:
% ai         [1x1]   initial semi-major axis            [km]
% ei         [1x1]   Initial eccentricity               [-]
% af         [1x1]   Final semi-major axis              [km]
% ef         [1x1]   Final eccentricity                 [-]
% type       [char]  maneuver type  
% mu         [1x1]   Gravitational parameter            [km^3/s^s]
%
% ---------------------------------------------------------------------
% Output arguments:
% DeltaV1       [1x1]   1st maneuver impulse            [km/s]
% DeltaV2       [1x1]   2nd maneuver impulse            [km/s]
% Deltat        [1x1]   maneuvre time                   [s]
%
% ----------------------------------------------------------------------

rp_i = a_i*(1-e_i);
ra_i = a_i*(1+e_i);
rp_f = a_f*(1-e_f);
ra_f = a_f*(1+e_f);

if type == "pa"
    
    rp_t = rp_i;
    ra_t = ra_f;

    a_t = (rp_t + ra_t)/2;
    e_t = (ra_t - rp_t)/(rp_t + ra_t);

    DeltaV1 = mu^(1/2)*((2/rp_t-1/a_t)^(1/2)-(2/rp_i-1/a_i)^(1/2));
    DeltaV2 = mu^(1/2)*((2/ra_f-1/a_f)^(1/2)-(2/ra_t-1/a_t)^(1/2));
    Deltat  = pi*((a_t^3)/mu)^(1/2);

elseif type == "ap"

    ra_t = ra_i;
    rp_t = rp_f;

    a_t = (rp_t + ra_t)/2;
    e_t = (ra_t - rp_t)/(rp_t + ra_t);

    DeltaV1 = mu^(1/2)*((2/ra_t-1/a_t)^(1/2)-(2/ra_i-1/a_i)^(1/2));
    DeltaV2 = mu^(1/2)*((2/rp_f-1/a_f)^(1/2)-(2/rp_t-1/a_t)^(1/2));
    Deltat  = pi*((a_t^3)/mu)^(1/2);

elseif type == "pp"

    rp_t = rp_i;
    ra_t = rp_f;

    a_t = (rp_t + ra_t)/2;
    e_t = (ra_t - rp_t)/(rp_t + ra_t);

    DeltaV1 = mu^(1/2)*((2/rp_t-1/a_t)^(1/2)-(2/rp_i-1/a_i)^(1/2));
    DeltaV2 = mu^(1/2)*((2/rp_f-1/a_f)^(1/2)-(2/ra_t-1/a_t)^(1/2));
    Deltat  = pi*((a_t^3)/mu)^(1/2);

elseif type == "aa"

    rp_t = ra_i;
    ra_t = ra_f;

    a_t = (rp_t + ra_t)/2;
    e_t = (ra_t - rp_t)/(rp_t + ra_t);

    DeltaV1 = mu^(1/2)*((2/rp_t-1/a_t)^(1/2)-(2/ra_i-1/a_i)^(1/2));
    DeltaV2 = mu^(1/2)*((2/ra_f-1/a_f)^(1/2)-(2/ra_t-1/a_t)^(1/2));
    Deltat  = pi*((a_t^3)/mu)^(1/2);

end








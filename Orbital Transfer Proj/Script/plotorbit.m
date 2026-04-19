function [X,Y,Z] = plotorbit(a, e, i, OM, om, th0, thf, dth, mu)

% 3D orbit plot
%
% [X,Y,Z] = plotorbit(a, e, i, OM, om, th0, thf, dth, mu)
% 
% -----------------------------------------------------------------------
% Input arguments:
% a         [1x1]   Semi-major axis             [km]
% e         [1x1]   Eccentricity                [-]
% i         [1x1]   Inclination                 [rad]
% OM        [1x1]   RAAN                        [rad]
% om        [1x1]   Pericenter anomaly          [rad]
% th0       [1x1]   Initial true anomaly        [rad]
% thf       [1x1]   Final true anomaly          [rad]
% dth       [1x1]   True anomaly step size      [rad]
% mu        [1x1]   Gravitational parameter     [km^3/s^s]
%
% ---------------------------------------------------------------------
% Output arguments:
% X         [dthx1] X coordinates of position vector    [km]
% Y         [dthx1] Y coordinates of position vector    [km]
% Z         [dthx1] Z coordinates of position vector    [km]
%
% ------------------------------------------------------------------------

% Se th0 o thf arrivano come vettori (es. da changePericenterArg), 
% prendiamo solo il primo valore.
th0 = th0(1); 
thf = thf(1);

% 2. GESTIONE SENSO DI PERCORRENZA:
% Se thf < th0 (es. da 350° a 10°), l'operatore ":" creerebbe un vettore vuoto.
% Aggiungiamo 2*pi per assicurarci che l'orbita giri nel verso giusto.
if thf < th0
    thf = thf + 2*pi;
end

% Creazione del vettore anomalie
th = th0 : dth : thf;

% 3. OTTIMIZZAZIONE (Pre-allocazione):
% Invece di r = [], pre-allochiamo per evitare rallentamenti.
n_punti = length(th);
r = zeros(n_punti, 3);

for n = 1:n_punti
    [rr, ~] = par2car(a, e, i, OM, om, th(n), mu);
    r(n, :) = rr'; % rr deve essere un vettore riga 1x3
end

% Output
X = r(:, 1);
Y = r(:, 2);
Z = r(:, 3);
end









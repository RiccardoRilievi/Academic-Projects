clear;
clc;


%% plot orbita finale

a_f = 14250.0;      % semi asse maggiore        [km]
e_f = 0.3205;       % eccentricità 
i_f = 0.9687;       % inclinazione              [rad]
OM_f = 0.3314;      % RAAN                      [rad]
om_f = 1.4320;      % anomalia del pericentro   [rad]
th_f = 2.9690;      % anomalia vera             [rad]

mu = 398600;

% analogamente all'orbita iniziale possiamo trovare il vettore velocità e
% posizione con la funziione par2car

[rr_f, vv_f]=par2car(a_f, e_f, i_f, OM_f, om_f, th_f, mu);

% altri parametri orbita finale

En_f = -mu / (2*a_f);           % Energia dell'orbita finale

rp_f = a_f * (1 - e_f);           % Raggio pericentro
ra_f = a_f * (1 + e_f);           % Raggio apocentro
p_f  = a_f * (1 - e_f^2);         % Semilatoretto
h_f  = (p_f * mu)^(1/2);        % h

vp_f = sqrt (mu / p_f) * (1 + e_f);       % velocità pericentro
va_f = sqrt (mu / p_f) * (1 - e_f);       % velocità apocentro

T_f = 2*pi * sqrt (a_f^3 / mu);         % Periodo

% e i vettori sul piano di riferimento equatoriale

th0 = th_f - 2*pi;
thf = th_f;
dth = 0.01;

[X_f,Y_f,Z_f] = plotorbit(a_f, e_f, i_f, OM_f, om_f, th0, thf, dth, mu);

% identifico pericentro

[Xp,Yp,Zp] = plotorbit(a_f, e_f, i_f, OM_f, om_f, 0, 0, dth, mu);

% identifico apocentro

[Xa,Ya,Za] = plotorbit(a_f, e_f, i_f, OM_f, om_f, pi, pi, dth, mu);

% plottiamo l'orbita iniziale

c1 = uisetcolor;
cp = uisetcolor;
ca = uisetcolor;

Terra_3D;

figure (1)
plot3(X_f,Y_f,Z_f, LineWidth=1, Color="#0072BD");
h = plot3 (nan, nan, nan, 'o', LineWidth=1, Color=c1);
set(h, 'XData', X_f(end), 'YData', Y_f(end), 'ZData', Z_f(end));
drawnow

plot3(Xp,Yp,Zp, "b", LineWidth=1);
h = plot3 (nan, nan, nan, 'o', LineWidth=1, Color=cp);
set(h, 'XData', Xp(end), 'YData', Yp(end), 'ZData', Zp(end));
drawnow

plot3(Xa,Ya,Za, "b", LineWidth=1);
h = plot3 (nan, nan, nan, 'o', LineWidth=1, Color=ca);
set(h, 'XData', Xa(end), 'YData', Ya(end), 'ZData', Za(end));
drawnow

grid on;
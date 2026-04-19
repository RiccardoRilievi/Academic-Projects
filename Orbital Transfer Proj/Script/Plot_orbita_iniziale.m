%% PLOT ORBITA 1

clear;
close all;
clc;

%% 1. CARATTERIZZAZIONE ORBITA INIZIALE

% Punto iniziale

rx = -7244.1296;    %
ry = -4751.0447;    %
rz = -465.5401;     % vettore posizione [km]

vx = 2.1140;        %
vy = -4.2990;       %
vz = -4.7890;       % vettore velocità [km/s]

rr = [rx; ry; rz];
vv = [vx; vy; vz];

mu = 398600;

[a, e, i, OM, om, th] = car2par(rr, vv, mu);  % trovati parametri orbitali dell'orbita di partenza

% altri parametri dell'orbita iniziale

En = -mu / (2*a);           % Energia dell'orbita iniziale

rp = a * (1 - e);           % Raggio pericentro
ra = a * (1 + e);           % Raggio apocentro
p  = a * (1 - e^2);         % Semilatoretto
h_i  = (p * mu)^(1/2);      % % h

vp = sqrt (mu / p) * (1 + e);       % velocità pericentro
va = sqrt (mu / p) * (1 - e);       % velocità apocentro

T = 2*pi * sqrt (a^3 / mu);         % Periodo

% Come prova dei risultati possiamo utilizzare la funzione inversa di
% car2par: par2car

[rr, vv]=par2car(a, e, i, OM, om, th, mu); % abbiamo appurato che i parametri trovati sono corretti

% plotorbit
%
% con la funzione plotorbit troviamo le componenti di ogni punto
% dell'orbita, per intervallo di tempo arbitrario, sul sistema di riferimento equatoriale
% 
% N.B 
% E' possibile utilizzare come th0 il th dell'orbita di partenza per
% valutare al meglio il punto di partenza, allo stesso modo è possibile utilizzare thf = th_f 
% per valutare il punto di arrivo 

th0 = th;
thf = 2*pi + th;
dth = 0.01;


[X,Y,Z] = plotorbit(a, e, i, OM, om, th0, thf, dth, mu);

% identifico pericentro

[Xp,Yp,Zp] = plotorbit(a, e, i, OM, om, 0, 0, dth, mu);

% identifico apocentro

[Xa,Ya,Za] = plotorbit(a, e, i, OM, om, pi, pi, dth, mu);

% plottiamo l'orbita iniziale

c1 = uisetcolor;
cp = uisetcolor;
ca = uisetcolor;

Terra_3D;

figure (1)
plot3(X,Y,Z, LineWidth=1, Color="#0072BD");
h = plot3 (nan, nan, nan, 'o', LineWidth=1, Color=c1);
set(h, 'XData', X(end), 'YData', Y(end), 'ZData', Z(end));
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
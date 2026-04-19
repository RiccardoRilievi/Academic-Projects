%% GRAFICO TRASFERIMENTO STANDARD

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

%% paramentri orbita finale

a_f = 14250.0;      % semi asse maggiore        [km]
e_f = 0.3205;       % eccentricità 
i_f = 0.9687;       % inclinazione              [rad]
OM_f = 0.3314;      % RAAN                      [rad]
om_f = 1.4320;      % anomalia del pericentro   [rad]
th_f = 2.9690;      % anomalia vera             [rad]
%--------------------------------------------------------

%%

th0 = th;
thf = 2*pi + th;
dth = 0.01;

[X,Y,Z] = plotorbit(a, e, i, OM, om, th0, thf, dth, mu);

% plottiamo l'orbita iniziale

Terra_3D;

figure (1)
plot3(X,Y,Z, "-.k");
h = plot3 (nan, nan, nan, 'og', LineWidth=1.5);
set(h, 'XData', X(end), 'YData', Y(end), 'ZData', Z(end));

drawnow
grid on;

% TRASFERIMENTO

% trasferimento bitangente (pa)

    % arrivo al pericentro 

    th1 = th;                                                                        
    th2 = 0;                                                                   
                                                                                 
    deltat_4_P = TOF (a, e, th1, th2, mu);

    [X,Y,Z] = plotorbit(a, e, i, OM, om, th, 2*pi, dth, mu);

figure (1)
plot3(X,Y,Z, "r", LineWidth=1);
h = plot3 (nan, nan, nan, 'or', LineWidth=1.5);
set(h, 'XData', X(end), 'YData', Y(end), 'ZData', Z(end));
drawnow

grid on;

type = "pa";

[DeltaV1_4_BT, DeltaV2_4_BT, Deltat_4_BT, a_t, e_t] = bitangetTransfer(a, e, a_f, e_f, type, mu);

DeltaV_4_BT = abs(DeltaV2_4_BT) + abs(DeltaV1_4_BT);

[X,Y,Z] = plotorbit(a_t, e_t, i, OM, om, -0.7, pi + 0.2, dth, mu);

figure (1)
plot3(X,Y,Z, "--k");
h = plot3 (nan, nan, nan);
set(h, 'XData', X(end), 'YData', Y(end), 'ZData', Z(end));
drawnow

[X,Y,Z] = plotorbit(a_t, e_t, i, OM, om, 0, pi, dth, mu);

figure (1)
plot3(X,Y,Z, "r", LineWidth=1);
h = plot3 (nan, nan, nan, 'or', LineWidth=1.5);
set(h, 'XData', X(end), 'YData', Y(end), 'ZData', Z(end));
drawnow

[X,Y,Z] = plotorbit(a_f, e_f, i, OM, om, 0, 2*pi, dth, mu);

figure (1)
plot3(X,Y,Z, "-.k");
h = plot3 (nan, nan, nan);
set(h, 'XData', X(end), 'YData', Y(end), 'ZData', Z(end));
drawnow

% cambio piano

[DeltaV_4_CP, omf_4_CP, theta_4_CP] = changeOrbitalPlane(a_f, e_f, i, OM, om, i_f, OM_f, mu);

th1 = pi;                                                                        
th2 = theta_4_CP;                                                                   
                                                                                 
deltat_4_CP = TOF (a_f, e_f, th1, th2, mu);

[X,Y,Z] = plotorbit(a_f, e_f, i, OM, om, 0, 2*pi, dth, mu);

figure (1)
plot3(X,Y,Z, "-.k");
h = plot3 (nan, nan, nan);
set(h, 'XData', X(end), 'YData', Y(end), 'ZData', Z(end));
drawnow

[X,Y,Z] = plotorbit(a_f, e_f, i, OM, om, pi, theta_4_CP, dth, mu);

figure (1)
plot3(X,Y,Z, "r", LineWidth=1);
h = plot3 (nan, nan, nan, "or", LineWidth=1.5);
set(h, 'XData', X(end), 'YData', Y(end), 'ZData', Z(end));
drawnow

% cambio anomalia del pericentro

[DeltaV_4_CPA, thi_4_CPA, thf_4_CPA] = changePericenterArg (a_f, e_f, omf_4_CP, om_f, mu);

[X,Y,Z] = plotorbit(a_f, e_f, i_f, OM_f, omf_4_CP, theta_4_CP - 0.2 , thi_4_CPA(1) + 0.5, dth, mu);
figure (1)
plot3(X,Y,Z, "--k");
h = plot3 (nan, nan, nan);
set(h, 'XData', X(end), 'YData', Y(end), 'ZData', Z(end));
drawnow

[X,Y,Z] = plotorbit(a_f, e_f, i_f, OM_f, omf_4_CP, theta_4_CP, thi_4_CPA(1), dth, mu);
figure (1)
plot3(X,Y,Z, "r", LineWidth=1);
h = plot3 (nan, nan, nan, "or", LineWidth=1.5);
set(h, 'XData', X(end), 'YData', Y(end), 'ZData', Z(end));
drawnow

th1 = theta_4_CP;                                                                  
th2 = thi_4_CPA(1);                                                                    
th3 = thi_4_CPA(2);                                                                    

deltat_CPA_1 = TOF (a_f, e_f, th1, th2, mu);                                       
deltat_CPA_2 = TOF (a_f, e_f, th1, th3, mu);                                        
                                                                                 
if deltat_CPA_1 < deltat_CPA_2                                                   
                                                                                  
     deltat_CPA = deltat_CPA_1;
        
else

    deltat_CPA = deltat_CPA_2;

end

deltat_4_CPA = deltat_CPA;

% arrivo al punto finale

th1 = thf_4_CPA(1);
th2 = th_f;

deltat_4_f = TOF (a_f, e_f, th1, th2, mu);

[X,Y,Z] = plotorbit(a_f, e_f, i_f, OM_f, om_f, thf_4_CPA(1), th_f, dth, mu);

figure (1)
plot3(X,Y,Z, "-.k");
h = plot3 (nan, nan, nan);
set(h, 'XData', X(end), 'YData', Y(end), 'ZData', Z(end));
drawnow

[X,Y,Z] = plotorbit(a_f, e_f, i_f, OM_f, om_f, thf_4_CPA, th_f, dth, mu);

figure (1)
plot3(X,Y,Z, "r", LineWidth=1);
h = plot3 (nan, nan, nan, "oc", LineWidth=1.5);
set(h, 'XData', X(end), 'YData', Y(end), 'ZData', Z(end));
drawnow

% delta v e delta t totali

DeltaT_4 = deltat_4_CP + deltat_4_CPA +...
            deltat_4_P + Deltat_4_BT + deltat_4_f

DeltaV_4 = abs(DeltaV_4_CP) + abs(DeltaV_4_CPA) + abs(DeltaV_4_BT)


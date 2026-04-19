
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

% plottiamo l'orbita iniziale


Terra_3D;

figure (1)
plot3(X,Y,Z, "-.k");
h = plot3 (nan, nan, nan, 'og', LineWidth=1);

grid on;

set(h, 'XData', X(end), 'YData', Y(end), 'ZData', Z(end));
drawnow




%% 2. CARATTERIZZAZIONE ORBITA FINALE

a_f = 14250.0;      % semi asse maggiore        [km]
e_f = 0.3205;       % eccentricità 
i_f = 0.9687;       % inclinazione              [rad]
OM_f = 0.3314;      % RAAN                      [rad]
om_f = 1.4320;      % anomalia del pericentro   [rad]
th_f = 2.9690;      % anomalia vera             [rad]

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

% plottiamo l'orbita iniziale


%figure(1)
%plot3(X_f,Y_f,Z_f);
%h = plot3 (nan, nan, nan, 'or');
%step_animation = 10;

%for j = 1:step_animation:length(X)

 %   set(h, 'XData', X_f(j), 'YData', Y_f(j), 'ZData', Z_f(j));
  %  drawnow
%end

%% STRATEGIA TRASFERIMENTO ALTERNATIVA 1 CP + CPA + SECANTE

% cambio piano

[DeltaV_1_CP, omf_1_CP, theta_1_CP] = changeOrbitalPlane(a, e, i, OM, om, i_f, OM_f, mu);

th1 = th;                                                                        
th2 = theta_1_CP;                                                                   
                                                                                 
deltat_1_CP = TOF (a, e, th1, th2, mu);

[X_f,Y_f,Z_f] = plotorbit(a, e, i, OM, om, th, theta_1_CP, dth, mu);

figure(1)
plot3(X_f,Y_f,Z_f, "r", LineWidth=1);
h = plot3 (nan, nan, nan, 'or', LineWidth=1);

    set(h, 'XData', X_f(end), 'YData', Y_f(end), 'ZData', Z_f(end));
    drawnow


[X_f,Y_f,Z_f] = plotorbit(a, e, i_f, OM_f, omf_1_CP, 0, 2*pi, dth, mu);

figure(1)
plot3(X_f,Y_f,Z_f, "--k");
h = plot3 (nan, nan, nan);

    set(h, 'XData', X_f(end), 'YData', Y_f(end), 'ZData', Z_f(end));
    drawnow


% cambio anomalia del pericentro

[DeltaV_1_CPA, thi_1_CPA, thf_1_CPA] = changePericenterArg (a, e, omf_1_CP, om_f, mu);

[X_f,Y_f,Z_f] = plotorbit(a, e, i_f, OM_f, omf_1_CP, theta_1_CP, thi_1_CPA(1), dth, mu);

figure(1)
plot3(X_f,Y_f,Z_f, "r", LineWidth=1);
h = plot3 (nan, nan, nan, "or", LineWidth=1);

    set(h, 'XData', X_f(end), 'YData', Y_f(end), 'ZData', Z_f(end));
    drawnow
    

th1 = theta_1_CP;                                                                  
th2 = thi_1_CPA(1);                                                                    
th3 = thi_1_CPA(2);                                                                    

deltat_CPA_1 = TOF (a, e, th1, th2, mu);                                       
deltat_CPA_2 = TOF (a, e, th1, th3, mu);                                        
                                                                                 
if deltat_CPA_1 < deltat_CPA_2                                                   
                                                                                  
     deltat_CPA = deltat_CPA_1;
        
else

    deltat_CPA = deltat_CPA_2;

end

deltat_1_CPA = deltat_CPA;


% Secante

    % trovo a ed e dell' orbita di trasferimento

    r_cpa = p   / ( 1 + e * cos(thf_1_CPA(1)));
    r_f   = norm(rr_f);

    e_t   = ( r_f - r_cpa) / ( r_cpa - r_f * cos (th_f - thf_1_CPA(1)))
    a_t   = r_cpa / ( 1 - e_t )
    p_t   = a_t * ( 1 - e_t^2 );

[X_f,Y_f,Z_f] = plotorbit(a_f, e_f, i_f, OM_f, om_f, th_f, th_f + 2*pi, dth, mu);

figure(1)
plot3(X_f,Y_f,Z_f, "-.k");
h = plot3 (nan, nan, nan, "oc", LineWidth=1);

    set(h, 'XData', X_f(end), 'YData', Y_f(end), 'ZData', Z_f(end));
    drawnow    

    [X_f,Y_f,Z_f] = plotorbit(a_t, e_t, i_f, OM_f, om_f + thf_1_CPA(1), -0.8, 3.25, dth, mu);

figure(1)
plot3(X_f,Y_f,Z_f, "--k");
h = plot3 (nan, nan, nan);

    set(h, 'XData', X_f(end), 'YData', Y_f(end), 'ZData', Z_f(end));
    drawnow

[X_f,Y_f,Z_f] = plotorbit(a_t, e_t, i_f, OM_f, om_f + thf_1_CPA(1), 0, th_f-thf_1_CPA(1), dth, mu);

figure(1)
plot3(X_f,Y_f,Z_f, "r", LineWidth=1);
h = plot3 (nan, nan, nan);

    set(h, 'XData', X_f(end), 'YData', Y_f(end), 'ZData', Z_f(end));
    drawnow

    % calcolo i delta v

    % punto 1 (trasferimento secante)

    v_rad_i = sqrt( mu / p ) * e * sin(thf_1_CPA(1));
    v_tr_i  = sqrt( mu / p ) * (1 + e * cos(thf_1_CPA(1)));

    v_tr_t1  = sqrt( mu / p_t ) * (1 + e_t);

    deltav_rad_1 = v_rad_i;
    deltav_tr_1  = v_tr_t1 - v_tr_i;

    DeltaV1      = sqrt( deltav_tr_1^2 + deltav_rad_1^2);

    % punto 2 (trasferimento secante)

    v_rad_f = sqrt( mu / p_f ) * e_f * sin(th_f);
    v_tr_f  = sqrt( mu / p_f ) * (1 + e_f * cos(th_f));

    v_rad_t2 = sqrt( mu / p_t ) * e_t * sin(th_f - thf_1_CPA(1));
    v_tr_t2  = sqrt( mu / p_t ) * (1 + e_t * cos(th_f - thf_1_CPA(1)));

    deltav_rad_2 = v_rad_f - v_rad_t2;
    deltav_tr_2  = v_tr_f - v_tr_t2;

    DeltaV2      = sqrt( deltav_tr_2^2 + deltav_rad_2^2);

    DeltaV_sec = abs(DeltaV1) + abs(DeltaV2)

    % calcolo il tempo di trasferimento

    th1 = 0;
    th2 = th_f - thf_1_CPA(1);

    deltat_sec = TOF (a_t, e_t, th1, th2, mu)

    % calcolo delta v e deltat totatli

    DeltaT_ALT1 = deltat_1_CPA + deltat_1_CP + deltat_sec;
    DeltaV_ALT1 = abs(DeltaV_1_CPA) + abs(DeltaV_1_CP) + abs(DeltaV_sec);

  










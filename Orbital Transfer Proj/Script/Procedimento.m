%% PROCEDIMENTO ELABORATO


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
plot3(X,Y,Z,"k", LineWidth=1);
h = plot3 (nan, nan, nan, '*r');
step_animation = 10;
grid on;

for j = 1:step_animation:length(X)

    set(h, 'XData', X(j), 'YData', Y(j), 'ZData', Z(j));
    drawnow
end



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


Terra_3D;

plot3(X_f,Y_f,Z_f);
h = plot3 (nan, nan, nan, 'or');
step_animation = 10;

for j = 1:step_animation:length(X)

    set(h, 'XData', X_f(j), 'YData', Y_f(j), 'ZData', Z_f(j));
    drawnow
end

%% 2. STRATEGIA DI TRASFERIMENTO

% innanzitutto possiamo ricordare che per le orbite circolari se:
% 
% § r_f / r_i < 11.94 si utilizza trasferimento alla Hohmann (Bitangente)
% §           > 15.58 si utilizza il biellittico  
% § se il valore è compreso invece dipende tutto dal raggio di apocentro
%   dell'orbita doi trasferimento

% 2.1 Scelta del tipo di trasferimento

% Alcuni esempi di valutazione delle verie strategie: 
%
% 1. Strategia per DeltaV minimo 
%       Se a_f/a_i >> conveniente trasferimento biellittico;
%       Il cambio di piano è meno costoso se lontano dall'attrattore, 
%       per esempio è possibile effettuare il cambio di piano durante il
%       trasferimento biellittico;
%       Potrebbe essere conveniente cambiare forma e poi cambiare piano.
%
% 2. Strategia per deltat minimo
%       Pericentro-Apocentro non è necessariamente il più breve;
%       Cambiare sequenza della manovra può ridurre le attese (cambia i
%       punti di manovra)
%
% N.B. CAMBIO DI PIANO SEMPRE PRIMA DEL CAMBIO DELL'ANOMALIA DEL
%      PERICENTRO!!!

% 2.? 

% 2.? Attesa fino al punto di manovra e cambio di piano

% [DeltaV_cOP, om_f, theta] = changeOrbitalPlane(a, e, i, OM, om, i_f, OM_f, mu);  % Con la function cOP, troviamo: 
                                                                                 % l'impulso necessario alla manovra (DeltaV)
% th1 = th;                                                                        % l'anomalia del pericentro dopo il cambio di piano (om_f)
% th2 = theta;                                                                     % l'anomalia vera del punto di manovra
                                                                                 %
% deltat_cOP = TOF (a, e, th1, th2, mu);                                           % Con la function TOF
                                                                                 % tempo di atteso dal punto di partenza al punto di manovra

                                               
% 2.? Attesa fino al punto di manovra e cambio anomalia del pericentro

% [DeltaV_cPA, thi, thf] = changePericenterArg (a, e, om, om_f, mu);               % function cPA, troviamo:
                                                                                 % l'impulso necessario alla manovra (DeltaV)
% th1 = theta;                                                                     % l'anomalia vera del punto di manovra:
% th2 = thi(2);                                                                    % due valori poichè i punti di intersezione delle ellissi sono due
% th3 = thf(2);                                                                    %
                                                                                 % da questi due valori valutiamo poi il tempo di attesa per effettuare 
% deltat_cPA_1 = TOF (a, e, th1, th2, mu);                                         % la manovra: poichè il costo, in termini di velocità, del cambio di anomalia 
% deltat_cPA_2 = TOF (a, e, th1, th3, mu);                                         % del pericentro non dipende dal punto in cui viene effettuato, possiamo 
                                                                                 % scegliere l'anomalia vera del punto di manovra in cui il tempo di attesa
% if deltat_cPA_1 < deltat_cPA_2                                                   % è minore
                                                                                 % 
%     deltat_cPA = deltat_cPA_1;
        
% else

%     deltat_cPA = deltat_cPA_2;

% end

%% 3 STRATEGIA DI TRASFERIMENTO

%% 3.1 STRATEGIA 1 CP+CPA+BT(pa) 

% cambio piano

[DeltaV_1_CP, omf_1_CP, theta_1_CP] = changeOrbitalPlane(a, e, i, OM, om, i_f, OM_f, mu);

th1 = th;                                                                        
th2 = theta_1_CP;                                                                   
                                                                                 
deltat_1_CP = TOF (a, e, th1, th2, mu);

% cambio anomalia del pericentro

[DeltaV_1_CPA, thi_1_CPA, thf_1_CPA] = changePericenterArg (a, e, omf_1_CP, om_f, mu);


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

% trasferimento bitangente (pa)

    % arrivo al pericentro 

    th1 = thf_1_CPA(1);                                                                        
    th2 = 0;                                                                   
                                                                                 
    deltat_1_P = TOF (a, e, th1, th2, mu);

type = "pa";

[DeltaV1_1_BT, DeltaV2_1_BT, Deltat_1_BT] = bitangetTransfer(a, e, a_f, e_f, type, mu);

DeltaV_1_BT = abs(DeltaV2_1_BT) + abs(DeltaV1_1_BT);

% arrivo al punto finale

th1 = pi;
th2 = th_f;

deltat_1_f = TOF (a_f, e_f, th1, th2, mu);

% delta v e delta t totali

DeltaT_1 = deltat_1_CP + deltat_1_CPA +...
            deltat_1_P + Deltat_1_BT + deltat_1_f;

DeltaV_1 = abs(DeltaV_1_CP) + abs(DeltaV_1_CPA) + abs(DeltaV_1_BT);

% grafico

figure(3)
plot(DeltaT_1, DeltaV_1,"o", LineWidth=5.5, color="#0072BD");
hold on;
grid on;

%% 3.2 STRATEGIA 2 CP+CPA+BT(ap)

% cambio piano

[DeltaV_2_CP, omf_2_CP, theta_2_CP] = changeOrbitalPlane(a, e, i, OM, om, i_f, OM_f, mu);

th1 = th;                                                                        
th2 = theta_2_CP;                                                                   
                                                                                 
deltat_2_CP = TOF (a, e, th1, th2, mu);

% cambio anomalia del pericentro

[DeltaV_2_CPA, thi_2_CPA, thf_2_CPA] = changePericenterArg (a, e, omf_2_CP, om_f, mu);


th1 = theta_2_CP;                                                                  
th2 = thi_2_CPA(1);                                                                    
th3 = thi_2_CPA(2);                                                                    

deltat_CPA_1 = TOF (a, e, th1, th2, mu);                                       
deltat_CPA_2 = TOF (a, e, th1, th3, mu);                                        
                                                                                 
if deltat_CPA_1 < deltat_CPA_2                                                   
                                                                                  
     deltat_CPA = deltat_CPA_1;
        
else

    deltat_CPA = deltat_CPA_2;

end

deltat_2_CPA = deltat_CPA;


% trasferimento bitangente (ap)

    % arrivo all'apocentro

    th1 = thf_1_CPA(1);                                                                        
    th2 = pi;                                                                   
                                                                                 
    deltat_2_A = TOF (a, e, th1, th2, mu);

type = "ap";

[DeltaV1_2_BT, DeltaV2_2_BT, Deltat_2_BT] = bitangetTransfer(a, e, a_f, e_f, type, mu);

DeltaV_2_BT = DeltaV2_2_BT + DeltaV1_2_BT;

% arrivo al punto finale

th1 = 0;
th2 = th_f;

deltat_2_f = TOF (a_f, e_f, th1, th2, mu);
    
% delta v e delta t totali

DeltaT_2 = deltat_2_CP + deltat_2_CPA +...
            deltat_2_A + Deltat_2_BT + deltat_2_f;

DeltaV_2 = abs(DeltaV_2_CP) + abs(DeltaV_2_CPA) + abs(DeltaV_2_BT);

% grafico

figure(3)
plot(DeltaT_2, DeltaV_2,"o", LineWidth=5.5, color="#D95319");

%% 3.3 STRATEGIA 3 CP+CPA+BE

% cambio piano

[DeltaV_3_CP, omf_3_CP, theta_3_CP] = changeOrbitalPlane(a, e, i, OM, om, i_f, OM_f, mu);

th1 = th;                                                                        
th2 = theta_3_CP;                                                                   
                                                                                 
deltat_3_CP = TOF (a, e, th1, th2, mu);

% cambio anomalia del pericentro

[DeltaV_3_CPA, thi_3_CPA, thf_3_CPA] = changePericenterArg (a, e, omf_3_CP, om_f, mu);


th1 = theta_3_CP;                                                                  
th2 = thi_3_CPA(1);                                                                    
th3 = thi_3_CPA(2);                                                                    

deltat_CPA_1 = TOF (a, e, th1, th2, mu);                                       
deltat_CPA_2 = TOF (a, e, th1, th3, mu);                                        
                                                                                 
if deltat_CPA_1 < deltat_CPA_2                                                   
                                                                                  
     deltat_CPA = deltat_CPA_1;
        
else

    deltat_CPA = deltat_CPA_2;

end

deltat_3_CPA = deltat_CPA;

% trasferimento biellittico

    % arrivo al pericentro

    th1 = thf_1_CPA(1);                                                                        
    th2 = 2*pi;                                                                   
                                                                                 
    deltat_3_P = TOF (a, e, th1, th2, mu);


ra_t = 30000;
[DeltaV1_3_BE, DeltaV2_3_BE, DeltaV3_3_BE, Deltat1_3_BE, Deltat2_3_BE] = biellipticTransfer(a, e, a_f, e_f, ra_t, mu);

DeltaV_3_BE = abs(DeltaV3_3_BE) + abs(DeltaV2_3_BE) + abs(DeltaV1_3_BE);
DeltaT_3_BE = Deltat1_3_BE + Deltat2_3_BE;

% arrivo al punto finale

th1 = 0;
th2 = th_f;

deltat_3_f = TOF (a_f, e_f, th1, th2, mu);

% delta v e delta t totali

DeltaT_3 = deltat_3_CP + deltat_3_CPA +...
            deltat_3_P + DeltaT_3_BE + deltat_3_f;

DeltaV_3 = abs(DeltaV_3_CP) + abs(DeltaV_3_CPA) + abs(DeltaV_3_BE);

% grafico

figure(3)
plot(DeltaT_3, DeltaV_3,"o", LineWidth=5.5, color="#EDB120");


%% 3.4 STRATEGIA 4 BT(pa)+CP+CPA

% trasferimento bitangente (pa)

    % arrivo al pericentro 

    th1 = th;                                                                        
    th2 = 0;                                                                   
                                                                                 
    deltat_4_P = TOF (a, e, th1, th2, mu);

type = "pa";

[DeltaV1_4_BT, DeltaV2_4_BT, Deltat_4_BT] = bitangetTransfer(a, e, a_f, e_f, type, mu);

DeltaV_4_BT = abs(DeltaV2_4_BT) + abs(DeltaV1_4_BT);

% cambio piano

[DeltaV_4_CP, omf_4_CP, theta_4_CP] = changeOrbitalPlane(a_f, e_f, i, OM, om, i_f, OM_f, mu);

th1 = pi;                                                                        
th2 = theta_4_CP;                                                                   
                                                                                 
deltat_4_CP = TOF (a_f, e_f, th1, th2, mu);

% cambio anomalia del pericentro

[DeltaV_4_CPA, thi_4_CPA, thf_4_CPA] = changePericenterArg (a_f, e_f, omf_4_CP, om_f, mu);


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

% delta v e delta t totali

DeltaT_4 = deltat_4_CP + deltat_4_CPA +...
            deltat_4_P + Deltat_4_BT + deltat_4_f;

DeltaV_4 = abs(DeltaV_4_CP) + abs(DeltaV_4_CPA) + abs(DeltaV_4_BT);

% grafico

figure(3)
plot(DeltaT_4, DeltaV_4,"o", LineWidth=5.5, color="#7E2F8E");


%% 3.5 STRATEGIA 5 BE+CP+CPA

% trasferimento biellittico

    % arrivo al pericentro

    th1 = th;                                                                        
    th2 = 0;                                                                   
                                                                                 
    deltat_5_P = TOF (a, e, th1, th2, mu);


ra_t = 30000;
[DeltaV1_5_BE, DeltaV2_5_BE, DeltaV3_5_BE, Deltat1_5_BE, Deltat2_5_BE] = biellipticTransfer(a, e, a_f, e_f, ra_t, mu);

DeltaV_5_BE = abs(DeltaV3_5_BE) + abs(DeltaV2_5_BE) + abs(DeltaV1_5_BE);
DeltaT_5_BE = Deltat1_5_BE + Deltat2_5_BE;

% cambio piano

[DeltaV_5_CP, omf_5_CP, theta_5_CP] = changeOrbitalPlane(a_f, e_f, i, OM, om, i_f, OM_f, mu);

th1 = 0;                                                                        
th2 = theta_5_CP;                                                                   
                                                                                 
deltat_5_CP = TOF (a_f, e_f, th1, th2, mu);

% cambio anomalia del pericentro

[DeltaV_5_CPA, thi_5_CPA, thf_5_CPA] = changePericenterArg (a_f, e_f, omf_4_CP, om_f, mu);


th1 = theta_5_CP;                                                                 
th2 = thi_5_CPA(1);                                                                    
th3 = thi_5_CPA(2);                                                                    

deltat_CPA_1 = TOF (a_f, e_f, th1, th2, mu);                                       
deltat_CPA_2 = TOF (a_f, e_f, th1, th3, mu);                                        
                                                                                 
if deltat_CPA_1 < deltat_CPA_2                                                   
                                                                                  
     deltat_CPA = deltat_CPA_1;
        
else

    deltat_CPA = deltat_CPA_2;

end

deltat_5_CPA = deltat_CPA;

% arrivo al punto finale

th1 = thf_5_CPA(1);
th2 = th_f;

deltat_5_f = TOF (a_f, e_f, th1, th2, mu);

% delta v e delta t totali

DeltaT_5 = deltat_5_CP + deltat_5_CPA +...
            deltat_5_P + DeltaT_5_BE + deltat_5_f;

DeltaV_5 = abs(DeltaV_5_CP) + abs(DeltaV_5_CPA) + abs(DeltaV_5_BE);

% grafico

figure(3)
plot(DeltaT_5, DeltaV_5,"o", LineWidth=5.5, color="#77AC30");

%% 3.6 STRATEGIA 6 CP+BT(pa)+CPA

% cambio piano

[DeltaV_6_CP, omf_6_CP, theta_6_CP] = changeOrbitalPlane(a, e, i, OM, om, i_f, OM_f, mu);

th1 = th;                                                                        
th2 = theta_6_CP;                                                                   
                                                                                 
deltat_6_CP = TOF (a, e, th1, th2, mu);

% trasferimento bitangente (pa)

    % arrivo al pericentro 

    th1 = theta_6_CP;                                                                        
    th2 = 0;                                                                   
                                                                                 
    deltat_6_P = TOF (a, e, th1, th2, mu);

type = "pa";

[DeltaV1_6_BT, DeltaV2_6_BT, Deltat_6_BT] = bitangetTransfer(a, e, a_f, e_f, type, mu);

DeltaV_6_BT = abs(DeltaV2_6_BT) + abs(DeltaV1_6_BT);

% cambio anomalia del pericentro

[DeltaV_6_CPA, thi_6_CPA, thf_6_CPA] = changePericenterArg (a_f, e_f, omf_6_CP, om_f, mu);

th1 = pi;                                                                  
th2 = thi_6_CPA(1);                                                                    
th3 = thi_6_CPA(2);                                                                    

deltat_CPA_1 = TOF (a_f, e_f, th1, th2, mu);                                       
deltat_CPA_2 = TOF (a_f, e_f, th1, th3, mu);                                        
                                                                                 
if deltat_CPA_1 < deltat_CPA_2                                                   
                                                                                  
     deltat_CPA = deltat_CPA_1;
        
else

    deltat_CPA = deltat_CPA_2;

end

deltat_6_CPA = deltat_CPA;

% arrivo al punto finale

th1 = thf_6_CPA(1);
th2 = th_f;

deltat_6_f = TOF (a_f, e_f, th1, th2, mu);

% delta v e delta t totali

DeltaT_6 = deltat_6_CP + deltat_6_CPA +...
            deltat_6_P + Deltat_6_BT+ deltat_6_f

DeltaV_6 = abs(DeltaV_6_CP) + abs(DeltaV_6_CPA) + abs(DeltaV_6_BT)

% grafico

figure(3)
plot(DeltaT_6, DeltaV_6,"o", LineWidth=5.5, color="#4DBEEE");

%% 3.7 STRATEGIA 7 CP+BE+CPA

% cambio piano

[DeltaV_7_CP, omf_7_CP, theta_7_CP] = changeOrbitalPlane(a, e, i, OM, om, i_f, OM_f, mu);

th1 = th;                                                                        
th2 = theta_7_CP;                                                                   
                                                                                 
deltat_7_CP = TOF (a, e, th1, th2, mu);

% trasferimento biellittico

    % arrivo al pericentro

    th1 = theta_7_CP;                                                                        
    th2 = 0;                                                                   
                                                                                 
    deltat_7_P = TOF (a, e, th1, th2, mu);


ra_t = 30000;
[DeltaV1_7_BE, DeltaV2_7_BE, DeltaV3_7_BE, Deltat1_7_BE, Deltat2_7_BE] = biellipticTransfer(a, e, a_f, e_f, ra_t, mu);

DeltaV_7_BE = abs(DeltaV3_7_BE) + abs(DeltaV2_7_BE) + abs(DeltaV1_7_BE);
DeltaT_7_BE = Deltat1_7_BE + Deltat2_7_BE;

% cambio anomalia del pericentro

[DeltaV_7_CPA, thi_7_CPA, thf_7_CPA] = changePericenterArg (a_f, e_f, omf_7_CP, om_f, mu);


th1 = 0;                                                                  
th2 = thi_7_CPA(1);                                                                    
th3 = thi_7_CPA(2);                                                                    

deltat_CPA_1 = TOF (a_f, e_f, th1, th2, mu);                                       
deltat_CPA_2 = TOF (a_f, e_f, th1, th3, mu);                                        
                                                                                 
if deltat_CPA_1 < deltat_CPA_2                                                   
                                                                                  
     deltat_CPA = deltat_CPA_1;
        
else

    deltat_CPA = deltat_CPA_2;

end

deltat_7_CPA = deltat_CPA;

% arrivo al punto finale

th1 = thf_7_CPA(1);
th2 = th_f;

deltat_7_f = TOF (a_f, e_f, th1, th2, mu);

% delta v e delta t totali

DeltaT_7 = deltat_7_CP + deltat_7_CPA +...
            deltat_7_P + DeltaT_7_BE+ deltat_7_f;

DeltaV_7 = abs(DeltaV_7_CP) + abs(DeltaV_7_CPA) + abs(DeltaV_7_BE);

% grafico

figure(3)
plot(DeltaT_7, DeltaV_7,"o", LineWidth=4.5, color="#A2142F");




%% STRATEGIA TRASFERIMENTO ALTERNATIVA 1 CP + CPA + SECANTE

% cambio piano

[DeltaV_A1_CP, omf_A1_CP, theta_A1_CP] = changeOrbitalPlane(a, e, i, OM, om, i_f, OM_f, mu);

th1 = th;                                                                        
th2 = theta_1_CP;                                                                   
                                                                                 
deltat_A1_CP = TOF (a, e, th1, th2, mu);

% cambio anomalia del pericentro

[DeltaV_A1_CPA, thi_A1_CPA, thf_A1_CPA] = changePericenterArg (a, e, omf_1_CP, om_f, mu);


th1 = theta_A1_CP;                                                                  
th2 = thi_A1_CPA(1);                                                                    
th3 = thi_A1_CPA(2);                                                                    

deltat_CPA_1 = TOF (a, e, th1, th2, mu);                                       
deltat_CPA_2 = TOF (a, e, th1, th3, mu);                                        
                                                                                 
if deltat_CPA_1 < deltat_CPA_2                                                   
                                                                                  
     deltat_CPA = deltat_CPA_1;
        
else

    deltat_CPA = deltat_CPA_2;

end

deltat_A1_CPA = deltat_CPA;

% Secante

    % trovo a ed e dell' orbita di trasferimento

    r_cpa = p   / ( 1 + e * cos(thf_A1_CPA(1)));
    r_f   = norm(rr_f);

    e_t   = ( r_f - r_cpa) / ( r_cpa - r_f * cos (th_f - thf_A1_CPA(1)));
    a_t   = r_cpa / ( 1 - e_t );
    p_t   = a_t * ( 1 - e_t^2 );

    % calcolo i delta v

    % punto 1 (trasferimento secante)

    v_rad_i = sqrt( mu / p ) * e * sin(thf_A1_CPA(1));
    v_tr_i  = sqrt( mu / p ) * (1 + e * cos(thf_A1_CPA(1)));

    v_tr_t1  = sqrt( mu / p_t ) * (1 + e_t);

    deltav_rad_1 = v_rad_i;
    deltav_tr_1  = v_tr_t1 - v_tr_i;

    DeltaV1_A1   = sqrt( deltav_tr_1^2 + deltav_rad_1^2);

    % punto 2 (trasferimento secante)

    v_rad_f = sqrt( mu / p_f ) * e_f * sin(th_f);
    v_tr_f  = sqrt( mu / p_f ) * (1 + e_f * cos(th_f));

    v_rad_t2 = sqrt( mu / p_t ) * e_t * sin(th_f - thf_1_CPA(1));
    v_tr_t2  = sqrt( mu / p_t ) * (1 + e_t * cos(th_f - thf_1_CPA(1)));

    deltav_rad_2 = v_rad_f - v_rad_t2;
    deltav_tr_2  = v_tr_f - v_tr_t2;

    DeltaV2_A1   = sqrt( deltav_tr_2^2 + deltav_rad_2^2);

    DeltaV_sec_A1 = abs(DeltaV1_A1) + abs(DeltaV2_A1);

    % calcolo il tempo di trasferimento

    th1 = 0;
    th2 = th_f - thf_1_CPA(1);

    deltat_sec = TOF (a_t, e_t, th1, th2, mu);

    % calcolo delta v e deltat totatli

    DeltaT_A1 = deltat_1_CPA + deltat_1_CP + deltat_sec;
    DeltaV_A1 = abs(DeltaV_A1_CPA) + abs(DeltaV_A1_CP) + abs(DeltaV_sec_A1);

    %c = uisetcolor;

    %figure(3)
    %plot(DeltaT_A1, DeltaV_A1,"o", LineWidth=2, Color=c);

%% 3.8 STRATEGIA 4 BT(ap)+CP+CPA

% trasferimento bitangente (pa)

    % arrivo all'apocentro 

    th1 = th;                                                                        
    th2 = pi;                                                                   
                                                                                 
    deltat_8_A = TOF (a, e, th1, th2, mu);

type = "ap";

[DeltaV1_8_BT, DeltaV2_8_BT, Deltat_8_BT] = bitangetTransfer(a, e, a_f, e_f, type, mu);

DeltaV_8_BT = abs(DeltaV2_8_BT) + abs(DeltaV1_8_BT);

% cambio piano

[DeltaV_8_CP, omf_8_CP, theta_8_CP] = changeOrbitalPlane(a_f, e_f, i, OM, om, i_f, OM_f, mu);

th1 = 0;                                                                        
th2 = theta_8_CP;                                                                   
                                                                                 
deltat_8_CP = TOF (a_f, e_f, th1, th2, mu);

% cambio anomalia del pericentro

[DeltaV_8_CPA, thi_8_CPA, thf_8_CPA] = changePericenterArg (a_f, e_f, omf_8_CP, om_f, mu);


th1 = theta_8_CP;                                                                  
th2 = thi_8_CPA(1);                                                                    
th3 = thi_8_CPA(2);                                                                    

deltat_CPA_1 = TOF (a_f, e_f, th1, th2, mu);                                       
deltat_CPA_2 = TOF (a_f, e_f, th1, th3, mu);                                        
                                                                                 
if deltat_CPA_1 < deltat_CPA_2                                                   
                                                                                  
     deltat_CPA = deltat_CPA_1;
        
else

    deltat_CPA = deltat_CPA_2;

end

deltat_8_CPA = deltat_CPA;

% arrivo al punto finale

th1 = thf_8_CPA(1);
th2 = th_f;

deltat_8_f = TOF (a_f, e_f, th1, th2, mu);

% delta v e delta t totali

DeltaT_8 = deltat_8_CP + deltat_8_CPA +...
            deltat_8_A + Deltat_8_BT + deltat_8_f;

DeltaV_8 = abs(DeltaV_8_CP) + abs(DeltaV_8_CPA) + abs(DeltaV_8_BT);

% grafico

figure(3)
plot(DeltaT_8, DeltaV_8,"bo", LineWidth=5.5);

%% 3.9 STRATEGIA 6 CP+BT(pa)+CPA

% cambio piano

[DeltaV_9_CP, omf_9_CP, theta_9_CP] = changeOrbitalPlane(a, e, i, OM, om, i_f, OM_f, mu);

th1 = th;                                                                        
th2 = theta_9_CP;                                                                   
                                                                                 
deltat_9_CP = TOF (a, e, th1, th2, mu);

% trasferimento bitangente (pa)

    % arrivo al pericentro 

    th1 = theta_9_CP;                                                                        
    th2 = pi;                                                                   
                                                                                 
    deltat_9_A = TOF (a, e, th1, th2, mu);

type = "ap";

[DeltaV1_9_BT, DeltaV2_9_BT, Deltat_9_BT] = bitangetTransfer(a, e, a_f, e_f, type, mu);

DeltaV_9_BT = abs(DeltaV2_9_BT) + abs(DeltaV1_9_BT);

% cambio anomalia del pericentro

[DeltaV_9_CPA, thi_9_CPA, thf_9_CPA] = changePericenterArg (a_f, e_f, omf_6_CP, om_f, mu);

th1 = pi;                                                                  
th2 = thi_9_CPA(1);                                                                    
th3 = thi_9_CPA(2);                                                                    

deltat_CPA_1 = TOF (a_f, e_f, th1, th2, mu);                                       
deltat_CPA_2 = TOF (a_f, e_f, th1, th3, mu);                                        
                                                                                 
if deltat_CPA_1 < deltat_CPA_2                                                   
                                                                                  
     deltat_CPA = deltat_CPA_1;
        
else

    deltat_CPA = deltat_CPA_2;

end

deltat_9_CPA = deltat_CPA;

% arrivo al punto finale

th1 = thf_9_CPA(1);
th2 = th_f;

deltat_9_f = TOF (a_f, e_f, th1, th2, mu);

% delta v e delta t totali

DeltaT_9 = deltat_9_CP + deltat_9_CPA +...
            deltat_9_A + Deltat_9_BT+ deltat_9_f

DeltaV_9 = abs(DeltaV_9_CP) + abs(DeltaV_9_CPA) + abs(DeltaV_9_BT)

% grafico

figure(3)
plot(DeltaT_9, DeltaV_9,"Ko", LineWidth=5.5);


    xlabel("Time (s)")
    ylabel("Velocity (km/s)")
    axis([10000 60000 2 4.5])
    grid minor
    legend("CP+CPA+BT(pa)", "CP+CPA+BT(ap)", "CP+CPA+BE", "BT(pa)+CP+CPA", ...
        "BE+CP+CPA", "CP+BT(pa)+CPA", "CP+BE+CPA","BT(ap)-CP-CPA", "CP-BT(ap)-CPA")



































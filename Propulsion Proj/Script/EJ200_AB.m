% Ciclo Termodinamico Turbofan a Flussi Associati
% Motore Eurojet EJ 200 regime supersonico AFTERBURNER ON

clc;
close all;
clear;


% Pedici:
% 0: ingresso alla presa d'aria
% 1: onda d'urto
% 2: uscita diffusore = ingresso fan
% 3: uscita fan (ramo caldo) = ingresso IGV
% igv: uscita IGV = ingresso compressore 
% 4: uscita compressore = ingresso camera combustione
% 5: uscita camera combustione = ingresso HPT
% 6: uscita HPT = ingresso LPT
% 7: uscita LPT = ingresso mixing
% 8: uscita mixing = ingresso after-burner
% 9: uscita after-burner = ingresso nozzle
% 10: uscita nozzle

%% Dati:

% Condizioni di volo
z = 10000; % [m] = [35000 ft]
m_a = 77; % [kg/s] portata ingresso

% Aria
cp_a = 1004; % [J/kgK]
gamma_a = 1.4;
R_a = 287; % [J/kgK]
cv_a = cp_a - R_a; % [J/kgK]

% Gas combusti (burner)
cp_gc = 1155; % [J/kgK]
gamma_gc = 1.33;
R_gc = 286.6; % [J/kgK]
cv_gc = cp_gc - R_gc; % [J/kgK]

% Diffusore
eta_d = 0.95;
D_d = 0.74; % [m]

% Fan
beta_f = 4.2;
BPR = 0.4;
eta_f = 0.88;
etam_f = 0.93;

% Compressore
beta_c = 6.2;
eta_c = 0.88;
etam_c = 0.93;

% Camera di combustione
eta_b = 0.98;
pi_b = 0.94;
T5_tot = 1900; % [K]
deltaH = 42e6; % [J/kg]

% Turbina alta pressione
eta_hpt = 0.92;
etam_hpt = 0.92;

% Turbina bassa pressione
eta_lpt = 0.92;
etam_lpt = 0.92;

% After-burner
eta_ab = 0.98;
pi_ab = 0.93;
gamma_ab = 1.30;
cp_ab = 1243; % [J/kgK]
R_ab = 286.6; % [J/kgK]
m_f_ab = 182 / 60; % [kg/s]

% Nozzle
eta_n = 0.96;

%% Procedimento

[T_a, rho_a, P_a] = atmosphere_vec(z); % Condizioni atmosferiche statiche ingresso

%% 0) Condizioni all'ingresso


a_in = sqrt(gamma_a * R_a * T_a);
v_in = 4*m_a / (rho_a * pi * D_d^2);
M0 = v_in / a_in;
T0_tot = T_a * (1 + (gamma_a - 1) * M0^2 / 2);
P0_tot = P_a * (1 + (gamma_a - 1) * M0^2 / 2) ^ (gamma_a / (gamma_a - 1));

S0 = 0; % [J/(kgK) ]

%% 1) Condizioni dopo onda d'urto

M1 = sqrt((M0^2 + 2 / (gamma_a - 1)) / (2 * gamma_a * M0 ^ 2 / (gamma_a - 1) - 1));

P1_tot = P0_tot * (((gamma_a + 1) * M0^2 / 2) / (1 + (gamma_a - 1) * M0^2 / 2)) ^ (gamma_a / (gamma_a - 1)) / (2 * gamma_a * M0^2 / (gamma_a + 1) - (gamma_a - 1) / (gamma_a + 1)) ^ (1 / (gamma_a - 1));
T1_tot = T0_tot;
P1 = P_a * (2 * gamma_a * M0^2 / (gamma_a + 1) - (gamma_a - 1) / (gamma_a + 1));
T1 = T_a * (1 + (gamma_a - 1) * M0^2 / 2) * (2 * gamma_a * M0^2 / (gamma_a - 1) - 1) * 2 * (gamma_a - 1) / ((gamma_a + 1)^2 * M0^2);
a1 = sqrt(gamma_a * R_a * T1);
v1 = M1 * a1;

dS1 = cp_a * log(T1_tot / T0_tot) - R_a * log(P1_tot / P0_tot);

%% 2) All'uscita della presa d'aria (ingresso del fan):

T2_tot = T1_tot;
T2_tot_id = eta_d * (T2_tot - T1) + T1;

P2_tot = P1 * (T2_tot_id / T1) ^ (gamma_a / (gamma_a - 1));
pi_d = P2_tot / P1_tot;

dS2 = dS1 + cp_a * log(T2_tot / T1_tot) - R_a * log(P2_tot / P1_tot);

%% 3) All'uscita del fan:

P3_tot = beta_f * P2_tot;

T3_tot_id = T2_tot * (beta_f) ^ ((gamma_a - 1) / gamma_a);
T3_tot = T2_tot + (T3_tot_id - T2_tot) / eta_f;

dS3 = dS2 + cp_a * log(T3_tot / T2_tot) - R_a * log(P3_tot / P2_tot);

%% IGV) VARIABLE INLET GUIDE VANES
M_igv=0.55;
pi_igv=0.98;
T_tot_igv=T3_tot;
P_tot_igv=pi_igv*P3_tot;

P_igv=P_tot_igv/((1+(gamma_a-1)/2*M_igv^2)^(gamma_a/(gamma_a-1)));
T_igv=T_tot_igv/((1+(gamma_a-1)/2*(M_igv^2)));

dS_igv = dS3 + cp_a * log(T_tot_igv / T3_tot) - R_a * log(P_tot_igv / P3_tot);

rho_igv=P_igv/(R_a*T_igv);
v_igv=M_igv*sqrt(gamma_a*R_a*T_igv);

%% 4) All'uscita del compressore:

P4_tot = beta_c * P3_tot;

T4_tot_id = T3_tot * (P4_tot / P3_tot) ^ ((gamma_a - 1) / gamma_a);
T4_tot = T3_tot + (T4_tot_id - T3_tot) / eta_c;

dS4 = dS_igv + cp_a * log(T4_tot / T_tot_igv) - R_a * log(P4_tot / P_tot_igv);

%% 5) All'uscita della camera di combustione:

T5_tot = 1900;
f_b = (cp_gc * T5_tot - cp_a * T4_tot) / (deltaH * eta_b - cp_gc * T5_tot);

P5_tot = pi_b * P4_tot;

dS5 = dS4 + cp_gc * log(T5_tot / T4_tot) - R_gc * log(P5_tot / P4_tot);

%% 6) All'uscita della HPT:

T6_tot = T5_tot - cp_a * (T4_tot - T3_tot) / (etam_c * etam_hpt * cp_gc * (1 + f_b));
T6_tot_id = T5_tot + (T6_tot - T5_tot) / eta_hpt;

P6_tot = P5_tot * (T6_tot_id / T5_tot) ^ (gamma_gc / (gamma_gc - 1));

dS6 = dS5 + cp_gc * log(T6_tot / T5_tot) - R_gc * log(P6_tot / P5_tot);

%% 7) All'uscita della LPT:

P7_tot = P3_tot; 

T7_tot_id = T6_tot * (P3_tot / P6_tot) ^ ((gamma_gc - 1) / gamma_gc);
T7_tot = eta_lpt * (T7_tot_id - T6_tot) + T6_tot;

dS7 = dS6 + cp_a * log(T7_tot / T6_tot) - R_a * log(P7_tot / P6_tot);

BPR = (etam_lpt * etam_f * cp_gc * (1 + f_b) * (T6_tot - T7_tot)) / (cp_a * (T3_tot - T2_tot)) - 1;

m_ah = m_a / (1 + BPR);
m_ac = m_ah * BPR;
m_f = f_b * m_ah;

cp_mix = (m_ac * cp_a + (m_ah + m_f) * cp_gc) / (m_a + m_f);
cv_mix = (m_ac * cv_a + (m_ah + m_f) * cv_gc) / (m_a + m_f);
gamma_mix = cp_mix / cv_mix;
R_mix = cp_mix - cv_mix;

%% 8) All'uscita della camera di mixing:

T8_tot = (BPR * cp_a * T3_tot + (1 + f_b) * cp_gc * T7_tot) / ((1 + f_b + BPR) * cp_mix);

P8_tot = P7_tot;

dS8 = dS7 + cp_mix * log(T8_tot / T7_tot) - R_mix * log(P8_tot / P7_tot);

%% 9) All'uscita dell'after-burner

T9_tot = (m_f_ab * deltaH * eta_ab + (m_a + m_f) * cp_mix * T8_tot) / ((m_a + m_f + m_f_ab) * cp_ab);

P9_tot = pi_ab * P8_tot;

dS9 = dS8 + cp_ab * log(T9_tot / T9_tot) - R_ab * log(P9_tot / P8_tot);

%% 10) All'uscita del nozzle

P10 = P_a;

T10_id = T9_tot * (P10 / P9_tot) ^ ((gamma_ab - 1) / gamma_ab);
T10 = T9_tot - eta_n * (T9_tot - T10_id);

dS10 = dS9 + cp_ab * log(T10 / T9_tot) - R_ab * log(P10 / P9_tot);

%% Condizioni all'uscita:

a_out = sqrt(gamma_ab * R_ab * T10);
v_out = sqrt(2 * cp_ab * (T9_tot - T10));
M_out = v_out / a_out;

%% costrusico i vettori di T e P

T_tot_vec    = [T0_tot; T1_tot; T2_tot; T3_tot; T_tot_igv; T4_tot; T5_tot; T6_tot; T7_tot; T8_tot];
T_tot_id_vec = [T0_tot; T1_tot; T2_tot_id; T3_tot_id; T4_tot_id; T5_tot; T6_tot_id; T7_tot_id; T8_tot; T9_tot; T10];
P_tot_vec    = [P0_tot; P1_tot; P2_tot; P3_tot; P4_tot; P5_tot; P6_tot; P7_tot; P8_tot; P9_tot; P10];
S_vec        = [S0; dS1; dS2; dS3; dS_igv; dS4; dS5; dS6; dS7; dS8];

T_tot_vec_cold = [T0_tot; T1_tot; T2_tot; T3_tot; T8_tot]; 
S_vec_cold     = [S0; dS1; dS2; dS3; dS8];

T_tot_h    = [T0_tot; T1_tot; T2_tot; T3_tot; T_tot_igv; T4_tot; T5_tot; T6_tot; T7_tot; T8_tot; T9_tot; T10];
T_tot_c    = [T0_tot; T1_tot; T2_tot; T3_tot; T3_tot; T3_tot; T3_tot; T3_tot; T3_tot; T8_tot; T9_tot; T10]; 

P_tot_h    = [P0_tot; P1_tot; P2_tot; P3_tot; P_tot_igv; P4_tot; P5_tot; P6_tot; P7_tot; P8_tot; P9_tot; P10];
P_tot_c    = [P0_tot; P1_tot; P2_tot; P3_tot; P3_tot; ; P3_tot; P3_tot; P3_tot; P3_tot; P8_tot; P9_tot; P10]; 


%% PROVA GRAFICI (SPERIAMO)

figure (1)
semilogy(S_vec, T_tot_vec, "ro-", LineWidth=1.25)
grid on;
hold on;
semilogy(S_vec_cold, T_tot_vec_cold, "bo-", LineWidth=1.25)
legend("flusso aria calda", "flusso aria fredda")

ylim([0 2500])
xlabel('Variazione entropia [J/kgK]')
ylabel('Variazione temperatura [K]')

X     = ["AMBIENTE", "INLET", "FAN","IGV", "HPC", "B-CHAMBER", "HPT", "LPT", "MIXING", "AB", "IN NOZZLE", "OUT NOZZLE"];
    
    C     = categorical(X, X);

    figure
    
    yyaxis left
    plot(C, T_tot_h, '-', LineWidth=1.25)
    grid on;
    hold on;
    plot(C, T_tot_c, '--', LineWidth=1.25)
    ylim([0 2500])
    ylabel('Temperatura [K]')

    yyaxis right
    plot(C, P_tot_h, '-', LineWidth=1.25)
    grid on;
    hold on;
    plot(C, P_tot_c, '--', LineWidth=1.25)
    ylabel('Pressione [Pa]')
        
    
    legend("T ramo caldo", "T ramo freddo", "P ramo caldo", "P ramo freddo")
    

%% PRESTAZIONI

T = (m_a + m_f + m_f_ab) * v_out - m_a * v_in; 

Is = T / m_a; 
Is_g = Is / 9.81; 

TSFC = ((m_f + m_f_ab) / T) *3600; 

P_p = T * v_in;

P_d = 0.5 * (m_a + m_f + m_f_ab) * (v_out - v_in)^2;

P_j = P_p + P_d;

P_a = (m_f + m_f_ab) * deltaH;

eta_th = P_j / P_a;

eta_pr = P_p / P_j;

eta_glob = eta_th * eta_pr;









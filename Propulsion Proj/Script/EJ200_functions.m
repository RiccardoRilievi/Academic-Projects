% Ciclo Termodinamico Turbofan a Flussi Associati
% Motore Eurojet EJ 200 regime supersonico

clc;
close all;
clear;

% Pedici:
% 0: ingresso alla presa d'aria
% 1: onda d'urto
% 2: uscita diffusore = ingresso fan
% 3: uscita fan (ramo caldo) = ingresso compressore
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
M0 = 1.2; % Mach all'ingresso

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
T5_tot = 1900.77; % [K]
deltaH = 42e6; % [J/kg]

% Turbina alta pressione
eta_hpt = 0.92;
etam_hpt = 0.92;

% Turbina bassa pressione
eta_lpt = 0.92;
etam_lpt = 0.92;

% After-burner

zip_AB = "ON"; % ON/OFF

eta_ab = 0.98;
pi_ab = 0.93;
gamma_ab = 1.30;
cp_ab = 1243; % [J/kgK]
R_ab = 286.6; % [J/kgK]
m_f_ab = 182 / 60; % [kg/s]

% Nozzle

zip_N = "EO"; % EO/CONV
eta_n = 0.96;

% Vettori grandezze vuoti

v_in_vec        = [];
m_a_vec         = [];
M1_vec          = [];
BPR_vec         = [];
M_out_vec       = [];
v_out_vec       = [];
a_out_vec       = [];

% Matrici vuote

T_tot_c_mat       = [];
P_tot_h_mat       = [];
P_tot_c_mat       = [];

prestazioni_mat   = [];

[T_a, rho_a, P_a] = atmosphere_vec(5000); % Condizioni atmosferiche statiche ingresso

%% Procedimento

for j=1:length(z)

    for i = 1:length(M0)
    
[T_a, rho_a, P_a] = atmosphere_vec(z(j)); % Condizioni atmosferiche statiche ingresso

%% 0) Condizioni all'ingresso

a_in = sqrt(gamma_a * R_a * T_a);
v_in = M0(i) * a_in;
T0_tot = T_a * (1 + (gamma_a - 1) * M0(i)^2 / 2);
P0_tot = P_a * (1 + (gamma_a - 1) * M0(i)^2 / 2) ^ (gamma_a / (gamma_a - 1));
m_a = rho_a * pi * v_in * D_d ^ 2 / 4;

S0 = 0; % [J/(molK) ]

%% 1) Condizioni dopo onda d'urto

[M1, T1_tot, P1_tot, T1, P1] = shockwave(M0(i), T0_tot, P0_tot, gamma_a, P_a, T_a);

dS1 = S0 + cp_a * log(T1_tot / T0_tot) - R_a * log(P1_tot / P0_tot);

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

%% 4) All'uscita del compressore:

P4_tot = beta_c * P3_tot;

T4_tot_id = T3_tot * (P4_tot / P3_tot) ^ ((gamma_a - 1) / gamma_a);
T4_tot = T3_tot + (T4_tot_id - T3_tot) / eta_c;

dS4 = dS3 + cp_a * log(T4_tot / T3_tot) - R_a * log(P4_tot / P3_tot);

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


%% 9) AFTERBURNER(ON, OFF) + NOZZLE + PRESTAZIONI

[T9_tot, P9_tot, T10_id, T10, P10, a_out, v_out, M_out, P_merito, D_out] = afterburner_nozzle(P_a, T8_tot, P8_tot, deltaH, cp_mix, gamma_mix, R_mix, eta_n, zip_AB, zip_N, m_a, m_f, v_in);

%% compilazione vettori e matrici

T_tot_h    = [T0_tot, T1_tot, T2_tot, T3_tot, T4_tot, T5_tot, T6_tot, T7_tot, T8_tot, T9_tot, T10];
T_tot_h_id = [T0_tot, T1_tot, T2_tot_id, T3_tot_id, T4_tot_id, T5_tot, T6_tot_id, T7_tot_id, T8_tot, T9_tot, T10_id];
T_tot_c    = [T0_tot, T1_tot, T2_tot, T3_tot, T3_tot, T3_tot, T3_tot, T3_tot, T8_tot, T9_tot, T10];
    
P_tot_h = [P0_tot, P1_tot, P2_tot, P3_tot, P4_tot, P5_tot, P6_tot, P7_tot, P8_tot, P9_tot, P10];
P_tot_c = [P0_tot, P1_tot, P2_tot, P3_tot, P3_tot, P3_tot, P3_tot, P3_tot, P8_tot, P9_tot, P10];

S_vec_h = [S0; dS1; dS2; dS3; dS4; dS5; dS6; dS7; dS8];
S_vec_c = [S0; dS1; dS2; dS3; dS3; dS3; dS3; dS3; dS8];

% GRAFICI

%X     = ["AMBIENTE", "INLET", "FAN", "HPC", "B-CHAMBER", "HPT", "LPT", "MIXING", "AB", "IN NOZZLE", "OUT NOZZLE"];
    
    %C     = categorical(X, X);

    %figure
    
    %yyaxis left
    %plot(C, T_tot_h, '-')
    %grid on;
    %hold on;
    %plot(C, T_tot_c, '--')

    %yyaxis right
    %plot(C, P_tot_h, '-')
    %grid on;
    %hold on;
    %plot(C, P_tot_c, '--')

    %legend("T ramo caldo", "T ramo freddo", "P ramo caldo", "P ramo freddo")

    %figure
    %semilogy(S_vec_h, T_tot_h(1:length(S_vec_h)), '-o');
    %grid on;
    %hold on;
    %semilogy(S_vec_c, T_tot_c(1:length(S_vec_c)), '-o');
    %legend("ramo caldo", "ramo freddo")

% Vettori

v_in_vec        = [v_in_vec; v_in];
m_a_vec         = [m_a_vec;  m_a];
M1_vec          = [M1_vec;   M1];
BPR_vec         = [BPR_vec;  BPR];
M_out_vec       = [M_out_vec; M_out];
v_out_vec       = [v_out_vec; v_out];
a_out_vec       = [a_out_vec; a_out];

% Matrici 

T_tot_c_mat       = [T_tot_c_mat; T_tot_c];
P_tot_h_mat       = [P_tot_h_mat; P_tot_h];
P_tot_c_mat       = [P_tot_c_mat; P_tot_c];

T_TOT_H (i, :, j) = T_tot_h;
prestazioni_mat (i, :, :) = P_merito;


    end
end 







   %% GRAFICI


    %X     = ["AMBIENTE", "INLET", "FAN", "HPC", "B-CHAMBER", "HPT", "LPT", "MIXING", "AB", "IN NOZZLE", "OUT NOZZLE"];
    
    %C     = categorical(X, X);

    %figure
    
    %yyaxis left
    %plot(C, T_tot_h, '-')
    %grid on;
    %hold on;
    %plot(C, T_tot_c, '--')

    %yyaxis right
    %plot(C, P_tot_h, '-')
    %grid on;
    %hold on;
    %plot(C, P_tot_c, '--')

    %legend("T ramo caldo", "T ramo freddo", "P ramo caldo", "P ramo freddo")

    figure
    semilogy(S_vec_h, T_tot_h(1:length(S_vec_h)), '-o');
    grid on;
    hold on;
    semilogy(S_vec_c, T_tot_c(1:length(S_vec_c)), '-o');
    legend("ramo caldo", "ramo freddo")
   
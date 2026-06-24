%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSIGNMENT 2 - STRUCTURAL DYNAMICS AND AEROELASTICITY
% Academic Year: 2025/2026
% Professors:    Prof. Andrea Zanoni 
% 
% Student:       Riccardo Rilievi
% ID:            10766323
%
% --- EXECUTION INSTRUCTIONS --- 
% Software:           MATLAB 
% Version:            R2025b (or later compatible version)
% Required Toolboxes: Symbolic Math Toolbox, Control System Toolbox
% Description: Ritz-Galerkin implementation for aeroelastic analysis,
%              including frequency convergence, frequency response function (FRF),
%              flutter/divergence speed optimization, and gust response.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
close all; 

%% DATA
syms y real
rho = 1.225;
b = 14;               % Wing span [m] 
ba = 2.5;
c = 1;                % Chord [m] 
gamma = deg2rad(15);  % Strut inclination angle [rad]
d = 0.15;
 
e = c/4;              % Distance AC to EA [m] 
CL_alpha = 2 * pi;    % Lift curve slope [1/rad]
CL_beta  = pi;
Cmac_beta = -0.1;
 
y_B = b/2;            % Strut attachment spanwise position [m]
x_B = - 0.5 * c;
m  = 27;
MT = 100;
EJ       = 1e7;                           % Bending stiffness [Nm^2]
GJ       = 1e6;                           % Torsional stiffness [Nm^2] 
EA_strut = 1e8;                           % Strut axial stiffness [N]
k_s      = 2 * EA_strut * cos(gamma) / b; 
J        = 7.56;                          % Moment of inertia per unit span [kg m^2 / m]
L1 = 9; L2 = 12;
yT = L1 * b / 26;
zT = -1 - 0.5 * L2 / 26;

%% PARAMETER STRUCT
p.rho      = rho;
p.b        = b;
p.ba       = ba;
p.c        = c;
p.gamma    = gamma;
p.d        = d;
p.e        = e;
p.CL_alpha = CL_alpha;
p.CL_beta  = CL_beta;
p.Cmac_beta= Cmac_beta;
p.y_B      = y_B;
p.m        = m;
p.MT       = MT;
p.EJ       = EJ;
p.GJ       = GJ;
p.k_s      = k_s;
p.J        = J;
p.yT       = yT;

%% RITZ-GALERKIN METHOD
N_max = 15;   % max number of iteration to study convergence
toll  = 1e-3; % relative error

% Vector initialization
omega1_vec = nan(N_max,1);   % 1st bending frequency history [rad/s]
omega2_vec = nan(N_max,1);   % 1st torsion frequency history [rad/s]
eps1_all   = nan(N_max,1);
eps2_all   = nan(N_max,1);

fprintf('%-4s  %-14s  %-14s  %-12s  %-12s\n', ...
    'N', 'omega_1 [rad/s]', 'omega_2 [rad/s]', 'eps_1', 'eps_2');
fprintf('------------------------------------------------------------------\n');

for i = 1 : N_max
    N = i;
    % Shape function definition
    Nt = N;
    Nw = Nt + 1;
    
    % Bending                  % Torsion
    phi = sym(zeros(Nw,1));    psi = sym(zeros(Nt, 1));
    
    for j = 1 : Nw
        phi(j) = (y/b)^(j+1);
    end
    
    for j = 1 : Nt
        psi(j) = (y/b)^j;
    end
    
    % Derivatives
    d2phi = diff(phi, y, 2);
    dpsi  = diff(psi, y, 1);
    
    %% MATRICES DEFINITION
    phi_yT = double(subs(phi, y, yT));
    phi_B  = double(subs(phi, y, y_B));
    psi_B  = double(subs(psi, y, y_B));
    psi_yT = double(subs(psi, y, yT));
    
    % Mass matrix
    Mww = double(int(m * phi * phi', y, 0, b) + MT * phi_yT * phi_yT');
    Mwt = double(int( - m * d * phi * psi', y, 0, b));
    Mtw = double(int( - m * d * psi * phi', y, 0, b));
    Mtt = double(int((m * d^2 + J) * psi * psi', y, 0, b) + MT * zT^2 * psi_yT * psi_yT'); 
    
    M = [Mww Mwt; Mtw Mtt];
    
    % Structural Stiffness Matrix
    Ks_ww = double(int(EJ * d2phi * d2phi', y, 0, b)) + (phi_B * phi_B') * sin(gamma)^2 * k_s ; 
    Ks_wt = -phi_B * psi_B' * x_B * sin(gamma)^2 * k_s;
    Ks_tw = -psi_B * phi_B' * x_B * sin(gamma)^2 * k_s;
    Ks_tt = double(int(GJ * dpsi * dpsi' , y, 0, b)) + (psi_B * psi_B')* x_B^2 * sin(gamma)^2 * k_s ;
    
    Ks = [Ks_ww Ks_wt; Ks_tw Ks_tt];
    
    % Aerodynamics Stiffness Matrix (multiplied by q_inf before use)
    Ka_ww = zeros(Nw);
    Ka_wt = double(int(c * CL_alpha * phi * psi', y, 0, b));
    Ka_tw = zeros(Nt, Nw);
    Ka_tt = double(int(e * c * CL_alpha * psi * psi', y, 0, b));
    
    Ka = [Ka_ww Ka_wt; Ka_tw Ka_tt];
    
    % Aerodynamics damping matrix (multiplied by q_inf before use)
    Ca_ww =   double(int(c * CL_alpha * phi * phi', y, 0, b));
    Ca_wt =   double(int(c * CL_alpha * phi * psi' * e, y, 0, b));
    Ca_tw =   double(int(c * CL_alpha * psi * phi' * e, y, 0, b));
    Ca_tt =   double(int(c * CL_alpha * psi * psi' * e^2, y, 0, b));
    
    Ca = [Ca_ww Ca_wt; Ca_tw Ca_tt];
    
    % RHS (gusts) (multiplied by vg * q_inf / U_inf before use)
    RHS_vg_1 = double(int(c * CL_alpha * phi, y, 0, b));
    RHS_vg_2 = double(int(c * CL_alpha * psi * e, y, 0, b));
    
    RHS_vg = [RHS_vg_1; RHS_vg_2];
    
    % RHS (beta) (multiplied by q_inf before use)
    RHS_beta_1 =  double(int(c * CL_beta * phi, y, b - ba, b));
    RHS_beta_2 =  double(int(e * c * CL_beta * psi + c^2 * Cmac_beta * psi, y, b - ba, b));
    
    RHS_beta = [RHS_beta_1; RHS_beta_2];

    % Eigenvalue Problem (Pure Structural)
    [V, D] = eig(full(double(Ks)), full(double(M))); 
    omega2 = diag(D); 
    valid = omega2 > 0;
    omega2 = omega2(valid);
    [omega2, ~] = sort(omega2); 
    omega = sqrt(omega2);
    omega1_vec(N) = omega(1);
    omega2_vec(N) = omega(2);
    
    if N == 1
        fprintf('%-4d  %-14.6f  %-14.6f  %-12s  %-12s\n', ...
            N, omega(1), omega(2), '—', '—');
    else
        eps1 = abs(omega1_vec(N)-omega1_vec(N-1)) / omega1_vec(N);
        eps2 = abs(omega2_vec(N)-omega2_vec(N-1)) / omega2_vec(N);
        eps1_all(N) = eps1;
        eps2_all(N) = eps2;
        fprintf('%-4d  %-14.6f  %-14.6f  %-12.2e  %-12.2e\n', ...
            N, omega(1), omega(2), eps1, eps2);
        
        if eps1 < toll && eps2 < toll
            converged = true;
            N_conv  = N;
            Nw_conv = Nw;
            Nt_conv = Nt;
        
            M_conv        = M;
            Ks_conv       = Ks;
            Ka_conv       = Ka;
            Ca_conv       = Ca;
            RHS_beta_conv = RHS_beta;
        
            phi_conv = phi;
            psi_conv = psi;
        
            fprintf('\n  >> Convergence achieved at N = %d (Nw=%d, Nt=%d) <<\n\n', ...
                N, Nw, Nt);
            break;
        end
    end                                                     
end


figure()
semilogy(1:N_max, eps1_all, '-o', 'LineWidth', 2, 'MarkerSize', 6)
hold on; grid on;
semilogy(1:N_max, eps2_all, '-s', 'LineWidth', 2, 'MarkerSize', 6)
yline(toll, 'r--', 'LineWidth', 2);
xlabel('Number of Modes N [-]', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Relative Frequency Error \epsilon [-]', 'FontSize', 12, 'FontWeight', 'bold')
title('Ritz-Galerkin Frequency Convergence History', 'FontSize', 13, 'FontWeight', 'bold')
legend('1^{st} Bending Mode (\epsilon_1)', '1^{st} Torsion Mode (\epsilon_2)', ...
       'Target Tolerance Threshold', 'FontSize', 11, 'Location', 'northeast')
set(gca, 'FontSize', 11)


%% FREQUENCY RESPONSE FUNCTION
% Solving the system M*q'' + Ca*q' + (Ks - Ka)*q = RHS_beta * beta(t)
U_inf = 25; % Airspeed [m/s]
q_inf = 0.5 * p.rho * U_inf^2;

% Compute the output at wing tip
xP = -0.5 * c;
phi_b  = double(subs(phi_conv, y, b));           
psi_b  = double(subs(psi_conv, y, b));           
output  = [phi_b ; -xP * psi_b];             

beta_0 = 1;
Omega_vec = logspace(-1, 3, 1000); 
H = zeros(1, length(Omega_vec));

for k = 1 : length(Omega_vec)
    Omega = Omega_vec(k);
    A     = - Omega^2 * M_conv + 1j * Omega * q_inf / U_inf * Ca_conv + (Ks_conv - q_inf * Ka_conv);
    Q = A \ (q_inf * RHS_beta_conv);
    H(k) = output' * Q / beta_0;
end


figure()

% Magnitude Plot
subplot(2,1,1)
semilogx(Omega_vec, 20*log10(abs(H)), '-', 'LineWidth', 2)
grid on; hold on;
xl1 = xline(omega1_vec(N_conv), 'r--', 'LineWidth', 1.5);
xl2 = xline(omega2_vec(N_conv), 'm--', 'LineWidth', 1.5);
xlabel('\Omega [rad/s]', 'FontSize', 11, 'FontWeight', 'bold')
ylabel('|H(j\Omega)| [dB]', 'FontSize', 11, 'FontWeight', 'bold')
title('FRF Magnitude Response (\beta_0 = 1 rad)', 'FontSize', 12, 'FontWeight', 'bold')
legend([xl1, xl2], {'\omega_1 (Bending)', '\omega_2 (Torsion)'}, 'FontSize', 10, 'Location', 'best')
set(gca, 'FontSize', 10)

% Phase Plot
subplot(2,1,2)
semilogx(Omega_vec, unwrap(angle(H))*180/pi, '-', 'LineWidth', 2)
grid on; hold on;
xline(omega1_vec(N_conv), 'r--', 'LineWidth', 1.5);
xline(omega2_vec(N_conv), 'm--', 'LineWidth', 1.5);
xlabel('\Omega [rad/s]', 'FontSize', 11, 'FontWeight', 'bold')
ylabel('Phase [deg]', 'FontSize', 11, 'FontWeight', 'bold')
title('FRF Phase Response', 'FontSize', 12, 'FontWeight', 'bold')
set(gca, 'FontSize', 10)



%% FLUTTER ANALYSIS
xB_vec = linspace(-0.5 * c, 0.5 * c, 100);
U_vec= linspace(0, 300, 500);
damping = nan(size(U_vec));
U_flutter = nan(size(xB_vec));
U_div = nan(size(xB_vec));

[Ks0, Ks1, Ks2, Ka0, Ca0] = precompute_matrices(p, phi_conv, psi_conv, y, Nw_conv, Nt_conv);

for i = 1:length(xB_vec)
    xB = xB_vec(i);
    Ks = Ks0  +  xB * Ks1  +  xB^2 * Ks2;
    eig_div = eig(Ks, Ka0);
    
    % Keep only real positive eigenvalues
    eig_div = eig_div(abs(imag(eig_div)) < 1e-8);
    eig_div = real(eig_div);
    eig_div = eig_div(eig_div > 0);
    if ~isempty(eig_div)
        U_div(i) = sqrt(min(eig_div));
    end
    
    max_re_prev = -inf;   
    tol = 1e-10;
    for j = 1:length(U_vec)
        U  = U_vec(j);
        Ka = U^2 * Ka0;
        Ca = U   * Ca0;
        n = size(M_conv, 1);
        A = [ zeros(n),            eye(n);
             -(Ks-Ka),              -Ca ];
        B = [eye(n), zeros(n); zeros(n) M_conv];
        lambda = eig(A, B);
        max_re = max(real(lambda));
        
        if j > 1 && max_re_prev < -tol && max_re >= tol
            U_flutter(i) = interp1([max_re_prev, max_re], [U_vec(j-1), U], 0);
            break;   
        end
        max_re_prev = max_re;
    end
end

[U_min, idx_min] = min(U_flutter);
[U_max, idx_max] = max(U_flutter);

fprintf('--- Stability Analysis ---\n')
fprintf('V_flutter at initial x_B = %.3f --> V_f = %.2f m/s\n', xB_vec(1)/c, U_flutter(1));
fprintf('Optimal x_B minimizing V_flutter: x_B/c = %.3f --> V_f = %.2f m/s\n\n', xB_vec(idx_min)/c, U_min);


figure()
plot(xB_vec/c, U_flutter, '-', 'LineWidth', 2.5)
hold on; grid on
plot(xB_vec/c, U_div, '-', 'LineWidth', 2.5)
m1 = xline(xB_vec(idx_min)/c, 'r--', 'LineWidth', 1.5);
m2 = xline(xB_vec(idx_max)/c, 'm--', 'LineWidth', 1.5);
xlabel('Normalized Chordwise Strut Attachment Position x_B / c [-]', 'FontSize', 12, 'FontWeight', 'bold'); 
ylabel('Critical Airspeed V [m/s]', 'FontSize', 12, 'FontWeight', 'bold')
title('Critical Flutter and Divergence Speeds vs Strut Position', 'FontSize', 13, 'FontWeight', 'bold')
ylim([0, 100]);
legend({'Flutter Speed V_{flutter}', 'Divergence Speed V_{divergence}', ...
        ['Min V_f Location (x_B/c = ' num2cnv(xB_vec(idx_min)/c) ')'], ...
        ['Max V_f Location (x_B/c = ' num2cnv(xB_vec(idx_max)/c) ')']}, ...
        'Location', 'best', 'FontSize', 11)
set(gca, 'FontSize', 11)



%% GUST RESPONSE WITH DIRECT RECOVERY
U_inf = 0.75 * U_flutter(1); 
q_inf = 0.5 * p.rho * U_inf^2;
vg_amp = 1; % Step gust amplitude [m/s]

n  = size(M_conv,1);
M  = double(M_conv);
Ks = double(Ks_conv);
Ka = double(q_inf * Ka_conv);
Ca = double((q_inf / U_inf) * Ca_conv);
Fg = double(RHS_vg * (q_inf / U_inf));

d2phi0 = zeros(Nw_conv,1);
for i = 1:Nw_conv
    d2phi = diff(phi_conv(i), y, 2);
    d2phi0(i) = double(subs(d2phi, y, 0));
end

C_Mx = p.EJ * [d2phi0' zeros(1,Nt_conv)];
s = tf('s');
Z_s = s^2 * M + s * Ca + (Ks - Ka);
Q_over_vG = minreal(Z_s \ Fg);
G_Mx = C_Mx * Q_over_vG;

Mx_statico = vg_amp * dcgain(G_Mx);   
t_samp = 10; 
[Mx_transitorio, t_out] = step(vg_amp * G_Mx, t_samp);

figure()
plot(t_out, Mx_transitorio, '-', 'LineWidth', 2)
hold on; grid on;
yline(Mx_statico, 'r--', 'LineWidth', 2)
xlabel('Time t [s]', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Root Bending Moment M_x(0,t) [Nm]', 'FontSize', 12, 'FontWeight', 'bold')
xlim([0, 4.5])
title('Direct Recovery Transient Response to a Unit Step Gust', 'FontSize', 13, 'FontWeight', 'bold')
legend('Transient Response', 'Asymptotic Static Value M_{x,static}', 'FontSize', 11, 'Location', 'best')
set(gca, 'FontSize', 11)


%% VARIANCE OF BENDING WITH LYAPUNOV EQUATION 
sys_Mx = ss((G_Mx));
Sv = 1; % White noise gust intensity [m^2/s^3]

% Direct state-space variance calculation
sigma2_Mx = covar(sys_Mx, Sv);
sigma_Mx  = sqrt(sigma2_Mx);

% Lyapunov Verification
A_ly = sys_Mx.A;
B_ly = sys_Mx.B;
C_ly = sys_Mx.C;

P = lyap(A_ly, B_ly * Sv * B_ly');
sigma2_lyap = C_ly * P * C_ly';

fprintf('--- Moment Variance ---\n')
fprintf('Root Bending Moment Variance:  sigma^2 = %.4e Nm^2\n', sigma2_Mx)
fprintf('Root Bending Moment Std. Dev:  sigma   = %.4e Nm\n', sigma_Mx)
fprintf('Lyapunov Equation Validation:  sigma^2 = %.4e Nm^2\n', sigma2_lyap)

%% functions

function [Ks0, Ks1, Ks2, Ka0, Ca0] = precompute_matrices(p, phi, psi, y, Nw, Nt)

d2phi = diff(phi, y, 2);
dpsi  = diff(psi,  y, 1);

phi_B = double(subs(phi, y, p.y_B));   % [Nw x 1]
psi_B = double(subs(psi, y, p.y_B));   % [Nt x 1]

% Structural Stiffness
Ks0_ww = double(int(p.EJ * d2phi * d2phi', y, 0, p.b)) + (phi_B * phi_B') * sin(p.gamma)^2 * p.k_s;
Ks0_tt = double(int(p.GJ * dpsi * dpsi', y, 0, p.b));
Ks0    = blkdiag(Ks0_ww, Ks0_tt);

Ks1_wt = - (phi_B * psi_B') * sin(p.gamma)^2 * p.k_s;
Ks1    = [zeros(Nw)      Ks1_wt;
    Ks1_wt'        zeros(Nt)];

Ks2_tt = (psi_B * psi_B') * sin(p.gamma)^2 * p.k_s;
Ks2    = blkdiag(zeros(Nw), Ks2_tt);

% Aerodynamic Stiffness
Ka0_wt = double(int(0.5*p.rho * p.c * p.CL_alpha * phi * psi',      y, 0, p.b));
Ka0_tt = double(int(0.5*p.rho * p.c * p.CL_alpha * psi * psi' * p.e, y, 0, p.b));

Ka0    = [zeros(Nw)      Ka0_wt; zeros(Nt, Nw)  Ka0_tt];

% Aerodynamic Damping
Ca0_ww = double(int(0.5*p.rho * p.c * p.CL_alpha        * phi * phi',        y, 0, p.b));
Ca0_wt = double(int(0.5*p.rho * p.c * p.CL_alpha        * phi * psi' * p.e,  y, 0, p.b));
Ca0_tt = double(int(0.5*p.rho * p.c * p.CL_alpha        * psi * psi' * p.e^2, y, 0, p.b));

Ca0 = [Ca0_ww   Ca0_wt; 
    Ca0_wt'  Ca0_tt];

end



function str = num2cnv(val)
% NUM2CNV Converts a numeric scalar into a formatted string.
%   str = num2cnv(val) takes a numeric value 'val' and returns a string 'str'
%   formatted to exactly 3 decimal places. This is used to maintain
%   constant notation inside plot legends and labels.
%
%
%   Inputs:
%       val - Numeric scalar value to be converted (e.g., 0.5 or -0.1534)
%
%   Outputs:
%       str - Formatted character array string (e.g., '0.500' or '-0.153')

% sprintf is a core MATLAB function that formats data into text.
% The format specifier '%.3f' forces a fixed-point notation with exactly 
% 3 digits after the decimal point, padding with zeros if necessary.
str = sprintf('%.3f', val);

end

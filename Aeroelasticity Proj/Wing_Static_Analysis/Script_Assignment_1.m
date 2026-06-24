%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSIGNMENT 1 - STRUCTURAL DYNAMICS AND AEROELASTICITY
% Academic Year: 2025/2026
% Professors: Prof. Giuseppe Quaranta, Prof. Andrea Zanoni 
% 
% Student: Riccardo Rilievi
% ID:      10766323
%
% --- EXECUTION INSTRUCTIONS --- 
% Software:          MATLAB 
% Version:           R2025b 
% Required Toolboxes: Symbolic Math Toolbox 
%
% Description: Numerical solution for the divergence analysis of a 
%              strut-braced wing using N Ritz shape functions.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;
close all; 

%% DATA

b = 14;               % Wing span [m] 
c = 1;                % Chord [m] 
EJ = 1e7;             % Bending stiffness [Nm^2]
GJ = 1e6;             % Torsional stiffness [Nm^2] 
EA_strut = 1e8;       % Strut axial stiffness [N] 
e = c/4;              % Distance AC to EA [m] 
CLa = 2*pi;           % Lift curve slope [1/rad]
gamma = deg2rad(15);  % Strut inclination angle [rad] 
y_B = b/2;            % Strut attachment spanwise position [m]
N = 5;               % Number of Ritz shape functions 

%% MATRICES DEFINITION
% Defining y as a symbolic variable 
syms y 

Kw = zeros(N); 
Ktheta = zeros(N);
Aw_theta = zeros(N); 
Atheta_theta = zeros(N);

for i = 1:N
    % Define shape functions for current index
    phi_i = (y/b)^(i+1);
    psi_i = (y/b)^i;
    
    for j = 1:N
        phi_j = (y/b)^(j+1);
        psi_j = (y/b)^j;
        
        % Stiffness Matrices 
        % Kw: Bending stiffness integral 
        Kw(i,j) = double(int(EJ * diff(phi_i, y, 2) * diff(phi_j, y, 2), y, 0, b));
        
        % Ktheta: Torsional stiffness integral
        Ktheta(i,j) = double(int(GJ * diff(psi_i, y, 1) * diff(psi_j, y, 1), y, 0, b));
        
        % Aerodynamic Matrices
        % Aw_theta: Coupling between torsion and bending displacement
        Aw_theta(i,j) = double(int(c * CLa * psi_j * phi_i, y, 0, b));
        
        % Atheta_theta: Aerodynamic moment due to torsion
        Atheta_theta(i,j) = double(int(e * c * CLa * psi_j * psi_i, y, 0, b));
    end
end

% Assembly of global wing matrices
K_wing = blkdiag(Kw, Ktheta);
A_aero = [zeros(N), Aw_theta; zeros(N), Atheta_theta];

%% 1.A DIVERGENCE DYNAMIC PRESSURE 

phi_B = zeros(N,1);
psi_B = zeros(N,1);
for i = 1:N
    phi_B(i) = (y_B/b)^(i+1);
    psi_B(i) = (y_B/b)^i;
end

% Strut equivalent vertical stiffness
L_strut = (y_B * tan(gamma)) / sin(gamma);
k_equiv = (EA_strut / L_strut) * (sin(gamma))^2;

% Chordwise offset study
offsets = linspace(-c/2, c/2, 100); 
q_div = zeros(size(offsets));

for k = 1:length(offsets)
    d = offsets(k);
    
    % Kinematic coupling vector at point B: h = [phi(yB); -d*psi(yB)] 
    h = [phi_B; -d*psi_B];
    
    % Total system stiffness including strut contribution
    K_total = K_wing + k_equiv * (h * h');
    
    % Solving the Generalized Eigenvalue Problem: [K]{x} = q[A]{x} 
    ev = eig(K_total, A_aero);
    
    % Filter for physical divergence pressure (lowest real positive eigenvalue)
    valid_ev = ev(isreal(ev) & ev > 0 & ~isinf(ev));
    
    if ~isempty(valid_ev)
        q_div(k) = min(valid_ev);
    else
        q_div(k) = NaN;
    end
end

% 1.A PLOT
figure('Color','w','Name','Divergence Analysis');
plot(offsets, q_div, 'b', 'LineWidth', 1.5);
grid on; hold on;

xline(-c/4, 'r--', 'AC');
xline(-c/2, 'k:', 'LE');
xline(c/2, 'k:', 'TE');
xline(0, 'k--', 'EA');

xlabel('Chordwise position of the strut to wing attachment wrt EA d [m]');
ylabel('Divergence dynamic pressure q_D [Pa]');
title('1.a: Divergence Dynamic Pressure VS Chordwise Position x_B');
legend('q_D','Aerodynamic Center', 'Leading Edge','Trailing Edge','Elastic Axis','Location','best');

%% 1.B CONTROL EFFECTIVENESS ANALYSIS

% 1.B.1 Aileron Geometric and Aerodynamic Data 
ba = 2.5;               % Aileron span [m]
CL_beta = pi;           % Aileron lift coefficient slope [1/rad]
CM_beta = -0.1;         % Aileron moment coefficient 
y_ail_start = b - ba;   % Aileron position: located at the wing tip 

% Aileron Load Vector {Q_beta}
% Calculated by integrating shape functions over the aileron span 
Q_w_beta = zeros(N,1);
Q_th_beta = zeros(N,1);

for i = 1:N
    
    phi_i = (y/b)^(i+1);
    psi_i = (y/b)^i;
    
    % Lift contribution due to aileron deflection
    Q_w_beta(i) = double(int(c * CL_beta * phi_i, y, y_ail_start, b));
    
    % Moment contribution around the Elastic Axis (EA)
    m_beta = (c^2 * CM_beta + e * c * CL_beta);
    Q_th_beta(i) = double(int(m_beta * psi_i, y, y_ail_start, b));
end
Q_beta_vec = [Q_w_beta; Q_th_beta];

% Rigid wing lift per unit dynamic pressure and unit deflection (L_R / (q * beta))
L_R_unit = c * CL_beta * ba;

% 1.B.3 Sensitivity Analysis Parameters & Dynamic Pressure Limits
gamma_vec = deg2rad(10:2:20);   % Strut inclination range [10:2:20] deg
d_fixed = -0.4;                 % Fixed strut offset [m] 
y_B = b/2;                      % Strut spanwise attachment point

% --- Dynamic Pressure Limit Calculation ---
% To ensure the plot range is physically relevant, we calculate the Divergence 
% Pressure (qD) specific to the chosen d_fixed and a reference gamma.

% 1. Calculate reference strut stiffness (using the median gamma for reference)
gamma_ref = mean(gamma_vec);
L_s_ref = (y_B * tan(gamma_ref)) / sin(gamma_ref);
k_eq_ref = (EA_strut / L_s_ref) * (sin(gamma_ref))^2;

% 2. Build the reference system matrix for the fixed offset
h_ref = [phi_B; -d_fixed * psi_B];
K_ref = K_wing + k_eq_ref * (h_ref * h_ref');

% 3. Find the divergence eigenvalues
ev_ref = eig(K_ref, A_aero);
valid_qD = ev_ref(isreal(ev_ref) & ev_ref > 0 & ~isinf(ev_ref));

if isempty(valid_qD)
    ev_r = eig(K_wing, A_aero);
    q_limit = 0.9 * min(ev_r(isreal(ev_r) & ev_r > 0));
    fprintf('Note: Configuration is inherently stable at d = %.2f. Using reference range.\n', d_fixed);
else
    q_limit = 0.95 * min(valid_qD);
end

q_range = linspace(1e-6, q_limit, 10000);
Ec_results = zeros(length(gamma_vec), length(q_range));

% Calculation

int_psi = zeros(1, N);
for i = 1:N
    psi_i = (y/b)^i;
    int_psi(i) = double(int(c * CLa * psi_i, y, 0, b));
end

for g = 1:length(gamma_vec)
    g_rad = gamma_vec(g);
    
    % Update strut stiffness for the current inclination gamma
    L_s = (y_B * tan(g_rad)) / sin(g_rad);
    k_eq = (EA_strut / L_s) * sin(g_rad)^2;
    
    % Kinematic coupling at point B with fixed offset
    h_b = [phi_B; -d_fixed * psi_B]; % Note the sign convention for offset d
    K_sys = K_wing + k_eq * (h_b * h_b');
    
    for jq = 1:length(q_range)
        q = q_range(jq);
        
        % Solve the Aeroelastic Equilibrium System: (K - qA)x = q{Q_beta} 
        
        state_sol = (K_sys - q * A_aero) \ (q * Q_beta_vec);
        
        % Extract torsional deformation coefficients (theta)
        theta_coeffs = state_sol(N+1:end);
        
        % Calculate extra lift due to elastic wing twist
        L_extra = q * (int_psi * theta_coeffs);
        
        % Control Effectiveness (E_c): Ratio of Total Lift to Rigid Lift 
        Ec_results(g, jq) = (q * L_R_unit + L_extra) / (q * L_R_unit);
    end
end

% 1.B PLOT
figure('Color','w','Name','Control Effectiveness - Reversal Focus');
hold on; grid on;
colors = lines(length(gamma_vec));

for g = 1:length(gamma_vec)
    plot(q_range, Ec_results(g,:), 'LineWidth', 2, 'Color', colors(g,:), ...
        'DisplayName', ['\gamma = ', num2str(round(rad2deg(gamma_vec(g)))), '^\circ']);
end


yline(0, 'k-', 'Reversal Condition','LineWidth', 1.5, 'HandleVisibility', 'off'); % Reversal line
yline(1, 'k:', 'Rigid Limit', 'HandleVisibility', 'off');

ylim([-0.5, 1.5]); 

xlabel('Dynamic pressure q [Pa]');
ylabel('Control Effectiveness E_C [-]');
title(['1.b: Aileron Effectiveness & Reversal (Offset d = ', num2str(d_fixed), ' m)']);
legend('Location', 'best');

%% 1.C COMPARISON OF DEFORMED SHAPES

% Test Parameters
beta0 = deg2rad(10);   % Aileron deflection [rad] (10 degrees)
d_ref = -0.4;           % Strut attachment [m]
gamma_ref = deg2rad(15); % Strut inclination [rad]

% Reference Dynamic Pressure Calculation
ev_clean = eig(K_wing, A_aero);
qD_clean = min(ev_clean(ev_clean > 0 & isreal(ev_clean)));
q_test = 0.4 * qD_clean; 

% Clean Wing (Without Strut)
% sol_clean contains modal participation factors [w_coeffs; theta_coeffs]
sol_clean = (K_wing - q_test * A_aero) \ (q_test * Q_beta_vec * beta0);

% Wing With Strut
L_s_ref = (y_B * tan(gamma_ref)) / sin(gamma_ref);
k_ref = (EA_strut / L_s_ref) * sin(gamma_ref)^2;
h_ref = [phi_B; -d_ref * psi_B];
K_strut_ref = k_ref * (h_ref * h_ref');
K_tot_strut = K_wing + K_strut_ref;

sol_strut = (K_tot_strut - q_test * A_aero) \ (q_test * Q_beta_vec * beta0);

% PHYSICAL DISPLACEMENTS
y_plot = linspace(0, b, 100);
w_clean = zeros(size(y_plot));     theta_clean = zeros(size(y_plot));
w_strut = zeros(size(y_plot));     theta_strut = zeros(size(y_plot));

for i = 1:N
    % Generate shape functions values along the span
    phi_y = (y_plot/b).^(i+1);
    psi_y = (y_plot/b).^i;
    
    % Superposition of modes for Clean Wing
    w_clean = w_clean + sol_clean(i) * phi_y;
    theta_clean = theta_clean + sol_clean(N+i) * psi_y;
    
    % Superposition of modes for Braced Wing
    w_strut = w_strut + sol_strut(i) * phi_y;
    theta_strut = theta_strut + sol_strut(N+i) * psi_y;
end

% 1.C PLOTS
figure('Color','w','Name','Task 1.c - Deformed Shapes Comparison');

subplot(2,1,1);
plot(y_plot, w_clean, 'r--', 'LineWidth', 2, 'DisplayName', 'No strut'); hold on;
plot(y_plot, w_strut, 'b', 'LineWidth', 2, 'DisplayName', 'With Strut');
ylabel('Bending displacement w [m]');
grid on;
title(['Static Aeroelastic Response at q = ', num2str(q_test, '%.1f'), ' Pa']);
legend('Location','northwest');

subplot(2,1,2);
plot(y_plot, rad2deg(theta_clean), 'r--', 'LineWidth', 2, 'DisplayName', 'No strut'); hold on;
plot(y_plot, rad2deg(theta_strut), 'b', 'LineWidth', 2, 'DisplayName', ' With Strut');
xlabel('Spanwise position y [m]');
ylabel('Torsion angle \theta [deg]');
grid on;
legend('Location','northwest');

%% 1.E CONVERGENCE ANALYSIS

N_max = N; 
N_vector = 1:N_max;
q_D_conv = nan(size(N_vector));
dqD_dN = zeros(size(N_vector)); 

gamma_conv = deg2rad(15);
d_conv = -0.4; 
L_s_conv = y_B / cos(gamma_conv);
k_eq_conv = (EA_strut / L_s_conv) * sin(gamma_conv)^2;


for n = N_vector
    
    idx = [1:n, (size(K_wing,1)/2 + 1):(size(K_wing,1)/2 + n)];
    
    K_sub = K_wing(idx, idx);
    A_sub = A_aero(idx, idx);
    
    h_sub = [phi_B(1:n); -d_conv * psi_B(1:n)];
    
    K_tot_sub = K_sub + k_eq_conv * (h_sub * h_sub');
    
    ev = eig(K_tot_sub, A_sub);
    
    valid = ev(isreal(ev) & ev > 1e-3 & ~isinf(ev));
    
    if ~isempty(valid)
        q_D_conv(n) = min(valid);
    end
end

figure('Color','w', 'Name', 'Task 1.e - Convergence Sensitivity');
plot(N_vector, q_D_conv, '-ok', 'LineWidth', 1.5, 'MarkerFaceColor', 'y', 'DisplayName', 'q_D Converging');
hold on; grid on; set(gca, 'GridLineStyle', ':');

yline(q_D_conv(end), '--r', sprintf('Asymptotic: %.1f Pa', q_D_conv(end)), ...
      'LineWidth', 1.5);

% Formatting
title('1.e: Convergence of the Ritz-Galerkin Solution', 'FontSize', 14);
ylabel('q_D [Pa]', 'FontSize', 12);
xlabel('Number of Shape Functions');
legend('Location', 'best');

%% 1. F STRUCTURAL REDESIGN 
fprintf('\n--- Starting Redesign & Sensitivity Analysis (Fixed d) ---\n');

% BASELINE CONFIGURATION 
d_fixed = -0.4; 
gamma_base = deg2rad(15);
L_s_base = y_B / cos(gamma_base);
h_base = [phi_B; -d_fixed * psi_B]; 
h_mat = h_base * h_base'; 

k_eq_base = (EA_strut / L_s_base) * (sin(gamma_base))^2;
ev_base = eig(K_wing + k_eq_base * h_mat, A_aero);
qD_baseline = min(ev_base(ev_base > 0 & isreal(ev_base)));
qD_target = 1.4 * qD_baseline;

fprintf('Baseline qD: %.2f Pa | Target qD (+40%%): %.2f Pa\n', qD_baseline, qD_target);

gamma_vec = deg2rad(10:1:40); 
fe_vec = 0.5:0.05:3.0;       
fw_vec = 1:0.1:1.5;    

mass_ratio_s_w = 0.05; 

Q_surf = zeros(length(fe_vec), length(gamma_vec));
best_mass_incr = inf;
found_sol = false;

for fw = fw_vec
    K_w_curr = K_wing * fw;
    m_w_incr_rel = (fw - 1); 
    
    for i = 1:length(fe_vec)
        fe = fe_vec(i);
        curr_EA = EA_strut * fe;
        
        for j = 1:length(gamma_vec)
            g = gamma_vec(j);
            L_s = y_B / cos(g);
            k_eq = (curr_EA / L_s) * (sin(g))^2;
            
            ev = eig(K_w_curr + k_eq * h_mat, A_aero);
            valid = ev(isreal(ev) & ev > 0);
            
            if isempty(valid)
                q_curr = NaN;
            else
                q_curr = min(valid);
            end
            
            if fw == 1.0
                Q_surf(i,j) = q_curr;
            end
            
            if ~isnan(q_curr) && q_curr >= qD_target
                m_s_incr_rel = (fe * L_s) / L_s_base - 1; 
                
                total_m_weighted = m_w_incr_rel + (mass_ratio_s_w * m_s_incr_rel);
                
                if total_m_weighted < best_mass_incr
                    best_mass_incr = total_m_weighted;
                    opt_res = struct('gamma', rad2deg(g), 'fe', fe, 'fw', fw, 'qD', q_curr);
                    found_sol = true;
                end
            end
        end
    end
end


% RESULTS
if found_sol
    fprintf('\n--- OPTIMAL REDESIGN FOUND ---\n');
    fprintf('Gamma: %.2f deg | fe: %.2f | fw: %.2f\n', opt_res.gamma, opt_res.fe, opt_res.fw);
    fprintf('Mass Increase: %.2f%%\n', best_mass_incr * 100);

else
    fprintf('\nNo solution found within the parameters.\n');
end
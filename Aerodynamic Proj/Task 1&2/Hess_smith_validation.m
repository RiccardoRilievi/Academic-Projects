%% Hess Smith Method
clc
clear 
addpath mat_functions

%% Input

U_inf = 1;      % Far-field velocity [m/s]
alpha = 1;      % Angle of attack
U_inf_x = U_inf * cos(deg2rad(alpha));
U_inf_y = U_inf * sin(deg2rad(alpha));
U_inf = [U_inf_x; U_inf_y];
Chord = 1;

LE_X_Position = 0;
LE_Y_Position = 0;

%% Create profile (with xfoil)

[NPanels] = extract_NPanels(['0008']);
[x,y]=createProfile(['0008'],NPanels,Chord);
geo.x=x;
geo.y=y;
figure(1)
plot(x,y,'-')
grid on
axis equal
hold on;
hold on

%% Create discretization & initialization
[centers,normals,tangent,extrema_1,extrema_2,alpha,lengths,L2G_TransfMatrix,G2L_TransfMatrix] = CreatePanels(geo);
        
NCols = sum(NPanels) + 1;
NRows = NCols;
A = zeros(NRows,NCols);     % system coefficients
B = zeros(NRows,1);         % known terms

%% Fill A
for i = 1:NPanels
    local_center = centers(i, :)';
    local_normal = normals(i, :)';
    for j = 1:NPanels
        local_extreme_1 = extrema_1(j, :)';
        local_extreme_2 = extrema_2(j, :)';
        local_L2G_TransfMatrix = squeeze(L2G_TransfMatrix(j, :, :));
        local_G2L_TransfMatrix = squeeze(G2L_TransfMatrix(j, :, :));
        A(i, j) = dot(uSource(local_center, local_extreme_1, local_extreme_2, local_L2G_TransfMatrix, local_G2L_TransfMatrix), local_normal);
        A(i, sum(NPanels)+1) = A(i, sum(NPanels)+1) + dot(uVortex(local_center, local_extreme_1, local_extreme_2, local_L2G_TransfMatrix, local_G2L_TransfMatrix), local_normal);
    end
end

%% Create a_v, c_s, and c_v vectors
first_centers = centers(1, :)';
first_tangent = tangent(1, :)';
last_centers = centers(end, :)';
last_tangent = tangent(end, :)';
last_a = 0;

for j = 1:NPanels
    local_extreme_1 = extrema_1(j, :)';
    local_extreme_2 = extrema_2(j, :)';
    local_L2G_TransfMatrix = squeeze(L2G_TransfMatrix(j, :, :));
    local_G2L_TransfMatrix = squeeze(G2L_TransfMatrix(j, :, :));
    a = dot(uSource(first_centers, local_extreme_1, local_extreme_2, local_L2G_TransfMatrix, local_G2L_TransfMatrix), first_tangent);
    last_a = last_a + dot(uVortex(first_centers, local_extreme_1, local_extreme_2, local_L2G_TransfMatrix, local_G2L_TransfMatrix), first_tangent);
    a = a + dot(uSource(last_centers, local_extreme_1, local_extreme_2, local_L2G_TransfMatrix, local_G2L_TransfMatrix), last_tangent);
    last_a = last_a + dot(uVortex(last_centers, local_extreme_1, local_extreme_2, local_L2G_TransfMatrix, local_G2L_TransfMatrix), last_tangent);
    A(sum(NPanels) + 1, j) = a;
end
A(sum(NPanels) + 1, sum(NPanels) + 1) = last_a;

%% Create B, the known terms in the system 
for j = 1:NPanels
    local_normal = normals(j, :)';
    B(j) = - dot(U_inf, local_normal);
end
first_tangent = tangent(1, :)';
last_tangent = tangent(end, :)';
B(sum(NPanels) + 1) = - dot(U_inf, (first_tangent + last_tangent));

%% Solve the linear system
solution = linsolve(A,B);

%% Compute velocity and cp
s_j = solution(1:NPanels);
v = solution(end);
U = zeros(NPanels,2);    
u_t = zeros(NPanels,1);   
u_n = zeros(NPanels,1);   
Cp = zeros(NPanels,1);
for i = 1:NPanels
    
    local_center = centers(i, :)';      
    t_i = tangent(i, :)';               
    n_i = normals(i, :)';               
    u_tot = U_inf;  
    for j = 1:NPanels
        local_extreme_1 = extrema_1(j, :)';
        local_extreme_2 = extrema_2(j, :)';
        local_L2G = squeeze(L2G_TransfMatrix(j, :, :));
        local_G2L = squeeze(G2L_TransfMatrix(j, :, :));
        u_s = uSource(local_center, local_extreme_1, local_extreme_2, local_L2G, local_G2L);
        u_v = uVortex(local_center, local_extreme_1, local_extreme_2, local_L2G, local_G2L);
        u_tot = u_tot + s_j(j) * u_s + v * u_v;   
        
    end
    u_t(i) = dot(u_tot, t_i);
    u_n(i) = dot(u_tot, n_i);
    Cp(i) = 1 - (u_t(i)/norm(U_inf))^2;
end
figure(4)
plot(centers(:, 1), Cp)
hold on;
grid on;
set(gca, 'YDir', 'reverse');  

%% Cl and Cm

Gamma = sum(solution(end) * lengths);
C_l = 2 * Gamma / norm(U_inf);
cl = [];
n_Uinf = [-U_inf_y; U_inf_x];

for i = 1:NPanels

cl_i = -dot(Cp(i) * lengths(i)' * normals(i, :), n_Uinf);
cl = [cl; cl_i];

end

C_l_alt = sum(cl);
cm = [];
centers_x_quart_chord = centers(:, 1) - 0.25;
centers_quart_chord = [centers_x_quart_chord centers(:, 2) zeros(length(centers),1)];
normals = [normals zeros(length(normals),1)];

for i = 1:NPanels

cm_i = dot(Cp(i) .* lengths(i) .* cross(centers_quart_chord(i,:), normals(i,:)), [0 0 -1]);
cm = [cm; cm_i];

end

C_M = -sum(cm);




%% matlab_results VS xfoil_results

cp_matlab = Cp;
x_matlab  = centers(:, 1);

xfoil_file = 'Cp_alpha1.dat';   % <-- CHANGE HERE!!!

data = readmatrix(xfoil_file);   
x_xfoil  = data(:,1);
cp_xfoil = data(:,2);


%% Plots comparison

figure; hold on; grid on;

% XFOIL
plot(x_xfoil, cp_xfoil, 'o-', 'LineWidth', 1.4, 'DisplayName','XFOIL');

% MATLAB
plot(x_matlab, cp_matlab, 'r-', 'LineWidth', 1.6, 'DisplayName','MATLAB');

set(gca,'YDir','reverse'); 
xlabel('x/c', 'FontSize', 20);
ylabel('C_p', 'FontSize', 20);
title('C_p - Hess Smith vs XFOIL', 'FontSize', 20);
legend('Location','best', 'FontSize', 20);
axis tight;


function [T_vec, rho_vec, p_vec] = atmosphere_vec(z_vec)

% ----------------------------------------------------------
%
% input
% 
% z_vec     = vector of altitude                      [ m ]
% ----------------------------------------------------------
%
% output
% T_vec     = vector temperature                      [ K ]
% rho_vec   = vector of density                       [ kg/m3 ]
% p_vec     = vector of pressure                      [ Pa ]
% ----------------------------------------------------------

R = 287.05;       % gas constant        [J/kgK]
g = 9.81;         %gravity acceleration [m/s^2]

p_vec   = [];
rho_vec = [];
T_vec   = [];

for i = 1 : length(z_vec)

if z_vec(i) <= 11e3

    lambda = -0.0065;   %[K/m]
    p_0 = 101325;       %[Pa]
    T_0 = 288.15;       %[K]
    rho_0 = 1.225;      %[kg/m^3]

    T = T_0 + lambda .* z_vec(i);
    rho = rho_0 * (T/T_0).^(- 1 - g / (R * lambda));
    p = p_0 * (T/T_0).^(-g / (R * lambda));

    else
    
    p_s = 22632;
    T_s = 216.55;
    rho_s = 0.363;
    h_s = 11e3;

    T = T_s;
    p = p_s * exp(-g .* (z_vec(i) - h_s) ./ (R * T_s));
    rho = rho_s * exp(-g .* (z_vec(i) - h_s) ./ (R * T_s));
end

    p_vec   = [p_vec; p];
    rho_vec = [rho_vec; rho];
    T_vec   = [T_vec; T];

end




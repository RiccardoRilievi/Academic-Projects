function [theta1, theta2]=tangentTransferOmega(a_i, e_i, a_f, e_f, x0)
%
% [theta1, theta2]=tangentTransferOmega(a_i, e_i, a_f, e_f, x0)
% -------------------------------------------------------------------------
% La funzione trova gli angoli tra il punto di tangenza di 2 orbite 
% rispetto al pericentro della prima e seconda orbita
% -------------------------------------------------------------------------
% Utilizzata per i trasferimenti tangenti previo 
% cambio argomento di pericentro
% -------------------------------------------------------------------------
% Input: 
% a_i        [1x1]   initial semi-major axis                [km]
% e_i        [1x1]   initial eccentricity                   [-]
% a_f        [1x1]   final semi-major axis                  [km]
% e_f        [1x1]   final eccentricity                     [-]
% x0         [1x1]   Initial guess for theta2               [rad]
% -------------------------------------------------------------------------
% Output:
% theta1     [1x1]   Angle between intersection and e1      [rad]
% theta2     [1x1]   Angle between intersection and e2      [rad]
% -------------------------------------------------------------------------
% Nota: per il calcolo è necessario che le 2 orbite abbiano 
%       le dimensioni adeguate per la condizione di tangenza
% -------------------------------------------------------------------------

p_i=a_i*(1-e_i^2);
p_f=a_f*(1-e_f^2);
r_a_i=p_i/(1-e_i);
r_p_f=p_f/(1+e_f);

% Condizione necessaria per la tangenza

if r_a_i<r_p_f
    error('Orbite non soddisfano la condizione di tangenza');
end

nmax=100;
toll=1e-16;

theta_1=@(theta_2) acos(p_i/(e_i*p_f)*(1+e_f*cos(theta_2))-1/e_i);
phi=@(theta_2) (e_i*sin(theta_1(theta_2)))/(1+e_i*cos(theta_1(theta_2)))-(e_f*sin(theta_2))/(1+e_f*cos(theta_2))+theta_2;

[succ, it] = ptofis(x0, phi, nmax, toll);

theta2=real(succ(end));
theta1=real(theta_1(theta2));

fprintf('Numero di iterazioni: %d\n', it);

end
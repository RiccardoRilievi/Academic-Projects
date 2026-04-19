function [T9_tot, P9_tot, T10_id, T10, P10, a_out, v_out, M_out, P_merito, D_out] = afterburner_nozzle(P_a, T8_tot, P8_tot, deltaH, cp_mix, gamma_mix, R_mix, eta_n, zip_AB, zip_N, m_a, m_f, v_in)

% After-burner

    eta_ab = 0.98;
    pi_ab = 0.93;
    gamma_ab = 1.30;
    cp_ab = 1243;       % [J/kgK]
    R_ab = 286.6;       % [J/kgK]
    m_f_ab = 182 / 60;  % [kg/s]

if zip_AB == "ON"

    T9_tot = (m_f_ab * deltaH * eta_ab + (m_a + m_f) * cp_mix * T8_tot) / ((m_a + m_f + m_f_ab) * cp_ab);

    P9_tot = pi_ab * P8_tot;

    % Uscita dal nozzle (EO/CONV)
    
    [T10, T10_id, P10, a_out, v_out, M_out] = nozzle (zip_N, P_a, T9_tot, P9_tot, gamma_ab, cp_ab, R_ab, eta_n);

    % prestazioni

    [T, Is, Is_g, TSFC, P_p, P_d, P_j, P_av, eta_th, eta_pr, eta_glob, D_out] = prestazioni(zip_AB, m_a, m_f, m_f_ab, T10, P10, P_a, R_ab, v_out, v_in, deltaH);

else

    P9_tot = P8_tot;
    T9_tot = T8_tot;

    % No AB

    % Uscita dal nozzle (EO/CONV)
    
    [T10, T10_id, P10, a_out, v_out, M_out] = nozzle (zip_N, P_a, T9_tot, P9_tot, gamma_mix, cp_mix, R_mix, eta_n);

    % Prestazioni

    [T, Is, Is_g, TSFC, P_p, P_d, P_j, P_av, eta_th, eta_pr, eta_glob, D_out] = prestazioni(zip_AB, m_a, m_f, m_f_ab, T10, P10, P_a, R_mix, v_out, v_in, deltaH);

end

P_merito = [T; Is; Is_g; TSFC; P_p; P_d; P_j; P_av; eta_th; eta_pr; eta_glob];

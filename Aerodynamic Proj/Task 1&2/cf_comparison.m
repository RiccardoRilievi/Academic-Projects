%% ====================== CF COMPARISON ======================= %%
clc; 
clear; 
close all;

%% ====== SETTINGS ======
folder = 'cf'; % folder with .dat files
AoA_target = 1;       % fixed angle of attack
N_fixed_compare = 8;  % N for 0008 vs 4417 comparison
%% ====== FILE READING ========

records = struct('airfoil',{}, 'N',{}, 'x',{}, 'Cf',{});

files = dir(fullfile(folder,'NACA*_cf_N*_A*.dat'));

% --- Read data and select upper surface ---
for k = 1:length(files)
    filename = files(k).name;
    filepath = fullfile(files(k).folder, filename);

    data = readmatrix(filepath);
    data = data(~all(isnan(data),2),:);

    if isempty(data), continue; end

    % Upper surface: take points up to maximum x
    [~, idx_maxx] = max(data(:,1));
    x  = data(1:idx_maxx,1);
    Cf = data(1:idx_maxx,2);

    % --- Limit to actual airfoil points (x <= 1) ---
    mask = x <= 1;
    x  = x(mask);
    Cf = Cf(mask);

    % Extract parameters from filename
    tok = regexp(filename,'NACA(\d+)_cf_N([\d\.]+)_A(\d+)','tokens');
    if isempty(tok), continue; end
    airfoil = tok{1}{1};
    N       = str2double(tok{1}{2});
    AoA     = str2double(tok{1}{3});

    if AoA ~= AoA_target, continue; end

    % Save record
    rec.airfoil = airfoil;
    rec.N       = N;
    rec.x       = x;
    rec.Cf      = Cf;
    records(end+1) = rec; %#ok<SAGROW>
end

if isempty(records)
    error('No files found for AoA = %d', AoA_target);
end

% Find unique airfoils
airfoils = unique(string({records.airfoil}));

%% ===== Cf(x) vs N =====

figure('Name','Cf vs N (upper surface)','Units','normalized','Position',[0.1 0.1 0.8 0.6]);
hold on; grid on;
legend_txt = {};

for a = 1:length(airfoils)
    af = airfoils{a};
    recs_af = records(strcmp(string({records.airfoil}), af));
    [~, order] = sort([recs_af.N]);
    recs_af = recs_af(order);

    for i = 1:length(recs_af)
        % Remove duplicates
        [x_unique, idx] = unique(recs_af(i).x);
        Cf_unique = recs_af(i).Cf(idx);

        if strcmp(af,'0008')
            style = '-';
        else
            style = '--';
        end

        plot(x_unique, Cf_unique, style, 'LineWidth',1.5)
        legend_txt{end+1} = sprintf('NACA %s, N=%.1f',af,recs_af(i).N);
    end
end

xlabel('x/c',  'FontSize', 20); 
ylabel('C_f (upper surface)', 'FontSize', 20);
title(sprintf('C_f distribution (upper surface) - AoA = %d°',AoA_target),  'FontSize', 20);
legend(legend_txt,'Location','best', 'FontSize', 10);

%% ===== Comparison 0008 vs 4417 for fixed N =====
figure('Name',sprintf('Cf vs x - comparison N=%g',N_fixed_compare),'Units','normalized','Position',[0.1 0.1 0.6 0.5]);
hold on; grid on;

for a = 1:length(airfoils)
    af = airfoils{a};
    recs_af = records(strcmp(string({records.airfoil}), af) & [records.N] == N_fixed_compare);

    if isempty(recs_af), continue; end

    % Remove duplicates
    [x_unique, idx] = unique(recs_af(1).x);
    Cf_unique = recs_af(1).Cf(idx);

    if strcmp(af,'0008')
        style = '-';
    else
        style = '--';
    end

    plot(x_unique, Cf_unique, style,'LineWidth',2)
end

xlabel('x/c'); 
ylabel('C_f (upper surface)');
title(sprintf('Comparison NACA 0008 vs 4417 - N=%g, AoA=%d°',N_fixed_compare,AoA_target));
legend('NACA 0008','NACA 4417','Location','best',  'FontSize', 15, 'FontWeight','bold');

%% ===== Separate plots for each airfoil =====
for a = 1:length(airfoils)
    af = airfoils{a};
    recs_af = records(strcmp(string({records.airfoil}), af));
    [~, order] = sort([recs_af.N]);
    recs_af = recs_af(order);

    figure('Name',sprintf('Cf vs x - NACA %s',af),'Units','normalized','Position',[0.1 0.1 0.6 0.5]);
    hold on; grid on;
    legend_txt = {};

    for i = 1:length(recs_af)
        % Remove duplicates
        [x_unique, idx] = unique(recs_af(i).x);
        Cf_unique = recs_af(i).Cf(idx);

        plot(x_unique, Cf_unique, 'LineWidth',2)
        legend_txt{end+1} = sprintf('N=%.1f', recs_af(i).N);
    end

    xlabel('x/c'); ylabel('C_f (upper surface)');
    title(sprintf('C_f distribution for NACA %s - AoA = %d°', af, AoA_target));
    legend(legend_txt,'Location','best');
end

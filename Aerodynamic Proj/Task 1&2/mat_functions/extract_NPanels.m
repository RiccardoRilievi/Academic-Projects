function [NPanels] = extract_NPanels(naca)

    % Read file .dat

        filename = strcat('NACA_', naca, '.dat');
        data = readmatrix(filename,'NumHeaderLines',1); % gestisce spazi multipli
        x = data(:,1);
        y = data(:,2);

        NPanels = length(x) - 1;

end
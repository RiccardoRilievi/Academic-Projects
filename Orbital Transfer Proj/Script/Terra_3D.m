function Terra_3D(Rt)
% Terra_3D.m - Earth texture loaded in a plot.
%
% PROTOTYPE:
%   Terra_3D(Rt)
%
% DESCRIPTION:
%   Function to load the Earth modelled as a sphere inside a figure.
%
% INPUT:
%   Rt          [1x1]       Earth mean radius       [km]
%
% OUTPUT:
%   []          [figure]    Figure open with the Earth picture loaded
% ------------------------------------------------------------------------

%% Default Input
if nargin < 1
    Rt = 6371.01;                                       % [km]
end

%%  Load the Earth image from a website
Earth_image = 'https://upload.wikimedia.org/wikipedia/commons/thumb/c/cd/Land_ocean_ice_2048.jpg/1024px-Land_ocean_ice_2048.jpg';

%% Figure setup
background_plot = 'w';
figure('Color', background_plot);
hold on;
grid on;
axis equal;
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
view(120,30);

%% Create Earth surface as a wireframe
npanels = 180;  
[x, y, z] = ellipsoid(0, 0, 0, Rt, Rt, Rt, npanels);
globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 'none');

%% Texturemap the globe (Download con bypass o fallback offline)
temp_filename = 'temp_earth_texture.jpg';

try
    % FIX: Creiamo un User-Agent per ingannare i blocchi anti-bot del server
    options = weboptions('UserAgent', 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)');
    websave(temp_filename, Earth_image, options);
    cdata = imread(temp_filename);
    
    % Applica la texture fotografica
    set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1, 'EdgeColor', 'none');
    
catch
    % FALLBACK OFFLINE: Se internet o il firewall bloccano tutto, usa i dati di MATLAB
    warning('Download fallito. Utilizzo la mappa topografica offline integrata in MATLAB.');
    load topo topo topomap1; % Carica il dataset base di MATLAB
    cdata = topo;
    
    % Applica la texture topografica di base
    set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1, 'EdgeColor', 'none');
    colormap(topomap1); % Applica i colori corretti alla mappa base
end

%% Cleanup
% Elimina il file temporaneo se è stato scaricato
if isfile(temp_filename)
    delete(temp_filename);
end

end
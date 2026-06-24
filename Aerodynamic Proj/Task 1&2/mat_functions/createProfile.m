function [x,y] = createProfile(naca,NPanels,Chord)
%   Detailed explanation goes here


        % Use xfoil to create the profile

        fileID = fopen('XFoilInput.txt','w');
        fprintf(fileID, ['naca ' ' '  naca, '\n\n']);
        fprintf(fileID,'pane\n\n');

        fprintf(fileID,'gdes\n');
        fprintf(fileID,'tgap 0 0 \n');
        
        fprintf(fileID,'exec \n\n\n');

        fprintf(fileID,'ppar\n');
        fprintf(fileID, ['n ' ' ' num2str(NPanels+1)  '\n\n\n']);

        filename = strcat('NACA_', naca, '.dat');

        fprintf(fileID, ['save ' ' ' filename '\n\n']);
        fprintf(fileID,'y\n\n');
        fprintf(fileID,'quit \n\n');
        fclose(fileID);

        Str2Exec = strcat("xfoil < XFoilInput.txt > /dev/null 2>&1");
%         Str2Exec = strcat("xfoil < XFoilInput.txt ");

        system(Str2Exec);

        % the import function was made with chatGPT!
        %Corpo = importXfoilProfile(filename);

         % Read file .dat
        data = readmatrix(filename,'NumHeaderLines',1); % gestisce spazi multipli
        x = data(:,1);
        y = data(:,2);

        % If necessary inverts vectors
        x = flipud(x) * Chord;
        y = flipud(y) * Chord;
end
        
       
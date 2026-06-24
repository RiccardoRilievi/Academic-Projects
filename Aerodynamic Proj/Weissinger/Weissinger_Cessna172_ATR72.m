close all
clear all %#ok<CLALL>
clc

Alpha = linspace(0, 3, 10); %[Deg]

% Cessna 172, Canadai CL-415 --------------------------------------------------------
PlaneName = ["Cessna172 Skyhawk","Canadair CL-415"];

% Wind Velocity
WindSpeed = [63.89, 92.5];

% Wing Parameters
RootChord       = [1.626; 3.53];
DihedralAngle   = [1.73; 0];
SweepAngle      = [2.08; 0];
TaperRatio      = [0.672; 1]; % RootChord / TipChord
Span            = [10.92; 28.60];

% Tail Parameters
RootChord_t     = [1.4; 2.64]; 
DihedralAngle_t = [0; 0];
SweepAngle_t    = [0; 0];
TaperRatio_t    = [0.57; 1];  % RootChord / TipChord
Span_t          = [3.4; 10.97];
Tail_x          = [4.5; 9.6];
Tail_y          = [0; 0];
Tail_z          = [-0.3; 1.8];      %minus?


% Preallocating For Speed
CL_To_Graph = zeros(length(Alpha), length(Span), 2);   %Number of different Alphae, Number of different Planes, Number of different Surfaces (Hardcoded to 2: Wing, Tail)
CD_To_Graph = zeros(length(Alpha), length(Span), 2);   %Number of different Alphae, Number of different Planes, Number of different Surfaces (Hardcoded to 2: Wing, Tail)
Total_CL_To_Graph = zeros(length(Alpha), length(Span));                     %iAlpha, iPlane
Total_CD_To_Graph = zeros(length(Alpha), length(Span));                     %iAlpha, iPlane

for iPlane  =1:length(Span)
    for iAlpha = 1:length(Alpha)

        %% Inputs

        U_Inf_Mag = WindSpeed(iPlane); 
        beta = 0;
        Alpha_U_inf = Alpha(iAlpha);
        U_Inf = [cosd(Alpha_U_inf)*cosd(beta) sind(beta) sind(Alpha_U_inf)] ./ norm([cosd(Alpha_U_inf) * cosd(beta) sind(beta) sind(Alpha_U_inf)]) .* U_Inf_Mag;
        rho = 1.225;
        
        config.NBodies = 2;
        
        config.RootChord     = [RootChord(iPlane),RootChord_t(iPlane)]; 
        config.DihedralAngle = [DihedralAngle(iPlane) ,DihedralAngle_t(iPlane)];     % [°]
        config.SweepAngle    = [SweepAngle(iPlane),SweepAngle_t(iPlane)];            % [°]
        config.TaperRatio    = [TaperRatio(iPlane), TaperRatio_t(iPlane)];           % TipChord/RootChord  
        config.Span          = [Span(iPlane),Span_t(iPlane)];
        config.LEPosition_X  = [0,Tail_x(iPlane)];
        config.LEPosition_Y  = [0,Tail_y(iPlane)];
        config.LEPosition_Z  = [0,Tail_z(iPlane)];
        
        config.RotationAngle_X = [0,0];   % ! Not Working !
        config.RotationAngle_Y = [0,0];   % ! Not Working !
        config.RotationAngle_Z = [0,0];   % ! Not Working !
        
        % Discretization options
        config.SemiSpanwiseDiscr = [20, 20];
        config.ChordwiseDiscr    = [20, 20];
        
        
        %% Preliminary computations
        
        % Computing the span
        config.SemiSpan = config.Span./2;
        % Computing the surface
        config.Surface = 2 * (config.SemiSpan .* config.RootChord .* ( 1 + config.TaperRatio ) ./ 2);
        config.SurfaceProjected = config.Surface .* cosd(config.DihedralAngle);
        % Computing the Tip chord
        config.TipChord = config.RootChord .* config.TaperRatio;
        
        % Compute MAC
        config.MAC = (2/3) .* config.RootChord .* ( (1 + config.TaperRatio + config.TaperRatio.^2)./(1 + config.TaperRatio));
        
        %% Create the geometry structure
        
        ControlPoints    = cell(config.NBodies, 1);
        InducedPoints    = cell(config.NBodies, 1);
        Normals          = cell(config.NBodies, 1);
        InfiniteVortices = cell(config.NBodies, 1);
        Vortices         = cell(config.NBodies, 1);
        internalMesh     = cell(config.NBodies, 1);
        WingExtremes     = cell(config.NBodies, 1);
        
        % Create L, D
        
        L_2Di       = cell(config.NBodies, 1);
        D_2Di       = cell(config.NBodies, 1);
        CL_2Di      = cell(config.NBodies, 1);
        CD_2Di      = cell(config.NBodies, 1);
        U_Section   = cell(config.NBodies, 1);
        L_3D        = cell(config.NBodies, 1);
        CL_3D       = cell(config.NBodies, 1);
        D_3D        = cell(config.NBodies, 1);
        CD_3D       = cell(config.NBodies, 1);
        
        for iBody = 1:config.NBodies
        
            [ControlPoints{iBody}, InducedPoints{iBody}, Normals{iBody}, InfiniteVortices{iBody}, Vortices{iBody}, internalMesh{iBody}, WingExtremes{iBody}] = createStructure(config, iBody);
        
        end
            
        
        %% Matrices initialization
        
        NPanelsTot = 0; 
        for iBody = 1:config.NBodies
            NPanelsTot = NPanelsTot+ 2* config.SemiSpanwiseDiscr(iBody) * config.ChordwiseDiscr(iBody)';
        end 
        
        matrixA = zeros(NPanelsTot, NPanelsTot);
        knownTerm = zeros(NPanelsTot, 1);
        
         %% Compute control point for each sector and sectors' Surface Area
        CPoint          = cell(config.NBodies, 1);
        SectorSurface   = cell(config.NBodies, 1);
        
        for iBody = 1:config.NBodies
            CPoint{iBody, 1} = zeros(config.SemiSpanwiseDiscr(iBody), 3);

            for iSector = 1:config.SemiSpanwiseDiscr(iBody) % spanwise only
                CPoint_LE = (internalMesh{iBody,1}{1,iSector}.LERoot + internalMesh{iBody,1}{1,iSector}.LEtip)/2;
                CPoint_TE = (internalMesh{iBody,1}{end,iSector}.TERoot + internalMesh{iBody,1}{end,iSector}.TEtip)/2;
            
                versor = (CPoint_LE - CPoint_TE) / norm(CPoint_LE - CPoint_TE);
            
                CPoint{iBody, 1}(iSector,:) = CPoint_TE + (3/4 * norm(CPoint_LE - CPoint_TE) ) * versor;

                BaseMinor   = norm(internalMesh{iBody,1}{1,iSector}.LEtip - internalMesh{iBody,1}{end,iSector}.TEtip);
                BaseMajor   = norm(internalMesh{iBody,1}{1,iSector}.LERoot - internalMesh{iBody,1}{end,iSector}.TERoot);
                Height      = internalMesh{iBody,1}{1,iSector}.LEtip(2) - internalMesh{iBody,1}{1,iSector}.LERoot(2);

                SectorSurface{iBody, 1}(1, iSector) = (BaseMajor + BaseMinor) .* Height / 2;
                
            end

            CPoint{iBody, 1} = [CPoint{iBody, 1}; flip([CPoint{iBody, 1}(:,1), -CPoint{iBody, 1}(:,2), CPoint{iBody, 1}(:,3)])];
            SectorSurface{iBody, 1} = [SectorSurface{iBody, 1}, flip(SectorSurface{iBody, 1})];
        end
        
        %% Construction of the matrix
        
        rowIndex = 0;
        
        for iBody = 1:config.NBodies
            
            % Cycle on all of its chordwise panels
            for ChordPanel_i = 1:config.ChordwiseDiscr(iBody)
                % Cycle on all of its spanwise panels
                for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iBody)
                    
                    % Update row index
                    rowIndex = rowIndex + 1;
           
                    columnIndex = 0;
                    
                    ControlPointHere = ControlPoints{iBody}{ChordPanel_i, SpanPanel_i}.Coords;
                    LocalNormal = Normals{iBody}{ChordPanel_i, SpanPanel_i}.Coords;
                    
                    
                    for jCorpo = 1:config.NBodies
                        
                        % Cycle on all of its chordwise panels
                        for ChordPanel_j = 1:config.ChordwiseDiscr(jCorpo)
                            % Cycle on all of its spanwise panels
                            for SpanPanel_j = 1:2*config.SemiSpanwiseDiscr(jCorpo)
                                
                                % Update column index
                                columnIndex = columnIndex + 1;
                                
                                % Compute the influence induced by first
                                % semi-infinite vortex
                                Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.toInfty;
                                Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.onWing;
                                U = vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);
                               
                                
                                % Compute the influence induced by finite vortex
                                Extreme_1 = Vortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root;
                                Extreme_2 = Vortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip;
                                U = U + vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);
                               
                                
                                % Compute the influence induced by second
                                % semi-infinite vortex
                                Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.onWing;
                                Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.toInfty;
                                U = U + vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);
                                
                                
                                matrixA(rowIndex, columnIndex) = dot(U, LocalNormal);
                                   
                            end
                        end
                    end
                end
            end
        end
        
        
        
        
        %% Costruzione del termine noto
        rowIndex = 0;
        for iBody = 1:config.NBodies
            
            % Cycle on all of its chordwise panels
            for ChordPanel_i = 1:config.ChordwiseDiscr(iBody)
                % Cycle on all of its spanwise panels
                for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iBody)
                    
                    % Update row index
                    rowIndex = rowIndex + 1;
          
                    LocalNormal = Normals{iBody}{ChordPanel_i, SpanPanel_i}.Coords;
                    
                    knownTerm(rowIndex) = -dot(U_Inf, LocalNormal);
                    
                end
            end
        end
        
        %% Solve the linear system
        
        Solution = linsolve(matrixA, knownTerm);
        
        Gamma = cell(config.NBodies, 1);
        
        rowIndex = 0;
        for iBody = 1:config.NBodies
            
            Gamma{iBody} = zeros( config.ChordwiseDiscr(iBody), config.SemiSpanwiseDiscr(iBody)*2 );
            
             % Cycle on all of its chordwise panels
            for ChordPanel_i = 1:config.ChordwiseDiscr(iBody)
                % Cycle on all of its spanwise panels
                for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iBody)
                    
                    % Update row index
                    rowIndex = rowIndex + 1;
                    
                    Gamma{iBody}(ChordPanel_i, SpanPanel_i) = Solution(rowIndex);
                end
            end
        end
        
        %% Induced velocity
        U_S = cell(config.NBodies, 1);
        U_Section_s = zeros(1,3);
        
        
        
        for iBody = 1:config.NBodies
                % Cycle on sectors
                U_S_Counter = 1;
                for Sectors_i = 1:2*config.SemiSpanwiseDiscr(iBody)
        
                    for jCorpo = 1:config.NBodies
        
                        % Cycle on all of its chordwise panels
                        for ChordPanel_j = 1:config.ChordwiseDiscr(jCorpo)
                            % Cycle on all of its spanwise panels
                            for SpanPanel_j = 1:2*config.SemiSpanwiseDiscr(jCorpo)
        
        
        
                                % Compute the influence induced by first
                                % semi-infinite vortex
                                Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.toInfty;
                                Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.onWing;
                                U_Section_s = U_Section_s + Gamma{jCorpo}(ChordPanel_j, SpanPanel_j)*vortexInfluence(CPoint{iBody,1}(Sectors_i,:), Extreme_1, Extreme_2);
        
        
        
                                % Compute the influence induced by second
                                % semi-infinite vortex
                                Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.onWing;
                                Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.toInfty;
                                U_Section_s = U_Section_s + Gamma{jCorpo}(ChordPanel_j, SpanPanel_j)*vortexInfluence(CPoint{iBody,1}(Sectors_i,:), Extreme_1, Extreme_2);
        
        
        
                            end
                        end
        
                    end

                   U_S{iBody, 1}(U_S_Counter, :) = U_Section_s;
                   U_S_Counter = U_S_Counter + 1;
                   U_Section_s = zeros(1,3);
        
                end
        end
        
        %% Compute the 2D and 3D Lift

        Delta_bi = cell(config.NBodies, 1); % Section Width
        
        for iBody = 1:config.NBodies
            Delta_bi{iBody, 1} = (config.SemiSpan(iBody)*2 / (config.SemiSpanwiseDiscr(iBody)*2)) * ones(config.SemiSpanwiseDiscr(iBody)*2, 1);

            L_2Di{iBody, 1} = zeros(1, config.SemiSpanwiseDiscr(iBody)*2);

            for iSector = 1:config.SemiSpanwiseDiscr(iBody)*2

                L_2Di{iBody}(1, iSector) = sum(rho * U_Inf_Mag * Gamma{iBody}(:, iSector) * cosd(config.DihedralAngle(iBody)));

            end
            CL_2Di{iBody} = L_2Di{iBody}(1, :)./( 0.5*rho*U_Inf_Mag^2 .* SectorSurface{iBody, 1}(:)');
        end
        

        for iBody = 1:config.NBodies

            L_3D{iBody, 1} = sum(Delta_bi{iBody, 1} .* L_2Di{iBody, 1}');
            CL_3D{iBody, 1} = L_3D{iBody, 1}./( 0.5 * rho * U_Inf_Mag^2 * config.Surface(1)); % Adimensionalized on the main Wing
        
        end
             

        %% Compute 2D and 3D induced drag
        
        alpha_i = cell(config.NBodies, 1);
        
        for iBody = 1:config.NBodies
            for iSector=1:config.SemiSpanwiseDiscr(iBody)*2
            
                    alphai = atan2( (dot(-U_S{iBody,1}(iSector,:), Normals{iBody,1}{1,iSector}.Coords)), U_Inf_Mag);
                    alpha_i{iBody,1} = [alpha_i{iBody,1}; alphai];
             
            end 
        
            D_2Di{iBody,1} = L_2Di{iBody,1}' .* sin(abs(alpha_i{iBody,1}));
            CD_2Di{iBody,1} = D_2Di{iBody,1}./( 0.5*rho*U_Inf_Mag^2 .* SectorSurface{iBody, 1}(:)); %always computed on the Main wing's surface
        
            D_3D{iBody,1} = sum(D_2Di{iBody,1}.*Delta_bi{iBody,1});
            CD_3D{iBody,1} = D_3D{iBody,1}./( 0.5*rho*U_Inf_Mag^2 *config.Surface(1));
        
        end 
        
        
        
        for iBody = 1:config.NBodies
            CL_To_Graph(iAlpha,iPlane,iBody) = CL_3D{iBody, 1};
            CD_To_Graph(iAlpha,iPlane,iBody) = CD_3D{iBody, 1};
            Total_CL_To_Graph(iAlpha, iPlane) = Total_CL_To_Graph(iAlpha, iPlane) + CL_To_Graph(iAlpha, iPlane, iBody); % This can be done as all Coefficients are normalized on the same Main Wing Area
            Total_CD_To_Graph(iAlpha, iPlane) = Total_CD_To_Graph(iAlpha, iPlane) + CD_To_Graph(iAlpha, iPlane, iBody); % This can be done as all Coefficients are normalized on the same Main Wing Area
        end 
        

    end

    %% -------------------------- Graphs -------------------------- Graphs -------------------------- Graphs --------------------------

    WingOrTail              = ["Wing", "Tail"];
    PlaneLineColor          = ["#1C1559", "#F5B427"]; % Cessna (Dark Blue), Canadair (Yellow)
    PlaneMarkerColor        = ["#FFFFFF", "#E62D02"]; % Cessna (White),     Canadair (Red)
    SurfaceMarkerType       = ["s", "o", "diamond"];  % Main Wing, Tail, Whole Aircraft
    SurfaceLineType         = ["-", "--", "-"];       % Main Wing, Tail, Whole Aircraft
    GraphsLineWidth         = 0.9;
    GraphsMarkerSize        = 5;
    TitleSize               = 20;
    AxisLabelSize           = 20;
    LegendSize              = 12;

    FigCounter = 1;

    % CL VS Alpha for each body (Wing, Tail)
    for iBody = 1:config.NBodies
        figure(FigCounter)
        plot(Alpha, CL_To_Graph(:,iPlane,iBody),'-o', 'DisplayName', sprintf('%s', PlaneName(iPlane)), 'Color', PlaneLineColor(iPlane), 'LineWidth', GraphsLineWidth, 'LineStyle', SurfaceLineType(iBody), 'MarkerFaceColor', PlaneMarkerColor(iPlane), 'Marker', SurfaceMarkerType(iBody), 'MarkerSize', GraphsMarkerSize)
        xlabel('\alpha [Deg]', 'FontSize', AxisLabelSize, 'FontWeight','bold')
        ylabel('C_{L, 3D} [-]', 'FontSize', AxisLabelSize, 'FontWeight','bold')
        hold on
        grid on

        if iBody == 1
            title('Wing C_{L, 3D} vs \alpha', '', 'FontSize', TitleSize)
        elseif iBody == 2
            title('Tail C_{L, 3D} vs \alpha', '', 'FontSize', TitleSize)
        end
        FigCounter = FigCounter + 1;
    end 

    % CL VS Alpha for whole Aircraft
    iBody = 3;
    figure(FigCounter)
    plot(Alpha, Total_CL_To_Graph(:,iPlane),'-o', 'DisplayName', sprintf('%s', PlaneName(iPlane)), 'Color', PlaneLineColor(iPlane), 'LineWidth', GraphsLineWidth, 'LineStyle', SurfaceLineType(iBody), 'MarkerFaceColor',PlaneMarkerColor(iPlane), 'Marker', SurfaceMarkerType(iBody), 'MarkerSize', GraphsMarkerSize)
    xlabel('\alpha [deg]', 'FontSize', AxisLabelSize)
    ylabel('C_{L, 3D} [-]', 'FontSize', AxisLabelSize)
    hold on
    grid on
    title('Aircraft C_{L, 3D} vs \alpha', '', 'FontSize', TitleSize)

    FigCounter = FigCounter + 1;

   
   % Polar for each body (Wing, Tail)
    for iBody = 1:config.NBodies
        figure(FigCounter)
        plot( CD_To_Graph(:,iPlane,iBody),  CL_To_Graph(:,iPlane,iBody),'-o','DisplayName',  sprintf('%s', PlaneName(iPlane)), 'Color', PlaneLineColor(iPlane), 'LineWidth', GraphsLineWidth, 'LineStyle', SurfaceLineType(iBody), 'MarkerFaceColor', PlaneMarkerColor(iPlane), 'Marker', SurfaceMarkerType(iBody), 'MarkerSize', GraphsMarkerSize)
        xlabel('C_{D, 3D}', 'FontSize', AxisLabelSize)
        ylabel('C_{L, 3D}', 'FontSize', AxisLabelSize)
        hold on
        grid on


        if iBody == 1
            title(sprintf('Wing Polar, \\alpha = %.1f', Alpha(iAlpha)), '', 'FontSize', TitleSize)
        elseif iBody == 2
            title(sprintf('Tail Polar, \\alpha = %.1f', Alpha(iAlpha)), '', 'FontSize', TitleSize)
        end

        FigCounter = FigCounter + 1;
    
    end

    % CD VS Alpha for whole Aircraft
    iBody = 3;
    figure(FigCounter)
    plot(Alpha, Total_CD_To_Graph(:,iPlane),'-o', 'DisplayName', sprintf('%s', PlaneName(iPlane)), 'Color', PlaneLineColor(iPlane), 'LineWidth', GraphsLineWidth, 'LineStyle', SurfaceLineType(iBody), 'MarkerFaceColor',PlaneMarkerColor(iPlane), 'Marker', SurfaceMarkerType(iBody), 'MarkerSize', GraphsMarkerSize)
    xlabel('\alpha [deg]', 'FontSize', AxisLabelSize)
    ylabel('C_{D, 3D} [-]', 'FontSize', AxisLabelSize)
    hold on
    grid on
    title('Aircraft C_{D, 3D} vs \alpha', '', 'FontSize', TitleSize)

    FigCounter = FigCounter + 1;


    % Polar for whole Aircraft
    iBody = 3;
    figure(FigCounter)
    plot(Total_CD_To_Graph(:,iPlane), Total_CL_To_Graph(:,iPlane),'-o', 'DisplayName', sprintf('%s', PlaneName(iPlane)), 'Color', PlaneLineColor(iPlane), 'LineWidth', GraphsLineWidth, 'LineStyle', SurfaceLineType(iBody), 'MarkerFaceColor', PlaneMarkerColor(iPlane), 'Marker', SurfaceMarkerType(iBody), 'MarkerSize', GraphsMarkerSize)
    xlabel('C_{D, 3D} [-]', 'FontSize', AxisLabelSize)
    ylabel('C_{L, 3D} [-]', 'FontSize', AxisLabelSize)
    hold on
    grid on
    title('Aircraft Polar', '', 'FontSize', TitleSize)

    FigCounter = FigCounter + 1;

    % L to D ratio VS Alpha (Wing, Tail)
    for iBody = 1:config.NBodies
        figure(FigCounter)
        plot(Alpha,  CL_To_Graph(:,iPlane,iBody)./CD_To_Graph(:,iPlane,iBody),'-o', 'DisplayName', sprintf('%s', PlaneName(iPlane)), 'Color', PlaneLineColor(iPlane), 'LineWidth', GraphsLineWidth, 'LineStyle', SurfaceLineType(iBody), 'MarkerFaceColor', PlaneMarkerColor(iPlane), 'Marker', SurfaceMarkerType(iBody), 'MarkerSize', GraphsMarkerSize)
        xlabel('\alpha [deg]', 'FontSize', AxisLabelSize)
        ylabel('C_{L}/C_{D} [-]', 'FontSize', AxisLabelSize)
        hold on
        grid on


        if iBody == 1
            title('Wing C_{L}/C_{D, induced}', '', 'FontSize', TitleSize)
        elseif iBody == 2
            title('Tail C_{L}/C_{D, induced}', '', 'FontSize', TitleSize)
        end
       FigCounter = FigCounter + 1;
    end 

    % CL Distribution (Wing, Tail)
    for iBody = 1:config.NBodies

        figure(FigCounter)
        plot(CPoint{iBody,1}(:,2),CL_2Di{iBody,1},'DisplayName', sprintf('%s %s', PlaneName(iPlane), WingOrTail(iBody)), 'Color', PlaneLineColor(iPlane), 'LineWidth', GraphsLineWidth, 'LineStyle', SurfaceLineType(iBody), 'MarkerFaceColor', PlaneMarkerColor(iPlane), 'Marker', SurfaceMarkerType(iBody), 'MarkerSize', GraphsMarkerSize)
        hold on
        grid on
        xlabel('Span [m]', 'FontSize', AxisLabelSize)
        ylabel('C_{L, 2D} [-]', 'FontSize', AxisLabelSize)
        title(sprintf('C_{L} Distributions, \\alpha = %.1f', Alpha(iAlpha)), '', 'FontSize', TitleSize)

    end 
    FigCounter = FigCounter + 1;

    % L Distribution (Wing, Tail)
    for iBody = 1:config.NBodies

        figure(FigCounter)
        plot(CPoint{iBody,1}(:,2), L_2Di{iBody,1},'DisplayName', sprintf('%s %s', PlaneName(iPlane), WingOrTail(iBody)), 'Color', PlaneLineColor(iPlane), 'LineWidth', GraphsLineWidth, 'LineStyle', SurfaceLineType(iBody), 'MarkerFaceColor', PlaneMarkerColor(iPlane), 'Marker', SurfaceMarkerType(iBody), 'MarkerSize', GraphsMarkerSize)
        hold on
        grid on
        xlabel('Span [m]', 'FontSize', AxisLabelSize)
        ylabel('L_{2D} [N]', 'FontSize', AxisLabelSize)
        title(sprintf('Lift Distributions, \\alpha = %.1f', Alpha(iAlpha)), '', 'FontSize', TitleSize)

    end 
    FigCounter = FigCounter + 1;

    % CD Distribution (Wing, Tail)
    for iBody = 1:config.NBodies

        figure(FigCounter)
        plot(CPoint{iBody,1}(:,2), CD_2Di{iBody,1},'DisplayName', sprintf('%s %s', PlaneName(iPlane), WingOrTail(iBody)), 'Color', PlaneLineColor(iPlane), 'LineWidth', GraphsLineWidth, 'LineStyle', SurfaceLineType(iBody), 'MarkerFaceColor', PlaneMarkerColor(iPlane), 'Marker', SurfaceMarkerType(iBody), 'MarkerSize', GraphsMarkerSize)
        hold on
        grid on
        xlabel('Span [m]', 'FontSize', AxisLabelSize)
        ylabel('C_{D, 2D} [-]', 'FontSize', AxisLabelSize)
        title(sprintf('C_{D} Distributions, \\alpha = %.1f', Alpha(iAlpha)), '', 'FontSize', TitleSize);

    end
    FigCounter = FigCounter + 1;

    % D Distribution (Wing, Tail)
    for iBody = 1:config.NBodies

        figure(FigCounter)
        plot(CPoint{iBody,1}(:,2), D_2Di{iBody,1},'DisplayName', sprintf('%s %s', PlaneName(iPlane), WingOrTail(iBody)), 'Color', PlaneLineColor(iPlane), 'LineWidth', GraphsLineWidth, 'LineStyle', SurfaceLineType(iBody), 'MarkerFaceColor', PlaneMarkerColor(iPlane), 'Marker', SurfaceMarkerType(iBody), 'MarkerSize', GraphsMarkerSize)
        hold on
        grid on
        xlabel('Span [m]', 'FontSize', AxisLabelSize)
        ylabel('D_{2D} [D]', 'FontSize', AxisLabelSize)
        title(sprintf('Drag Distributions, \\alpha = %.1f', Alpha(iAlpha)), '', 'FontSize', TitleSize);

    end
    FigCounter = FigCounter + 1;

    % Subplot CL and CD Distributions
    figure(FigCounter)
    ax1 = subplot(1, 2, 1);
    for iBody = 1:config.NBodies % CL

        plot(CPoint{iBody,1}(:,2),CL_2Di{iBody,1},'DisplayName', sprintf('%s %s', PlaneName(iPlane), WingOrTail(iBody)), 'Color', PlaneLineColor(iPlane), 'LineWidth', GraphsLineWidth, 'LineStyle', SurfaceLineType(iBody), 'MarkerFaceColor', PlaneMarkerColor(iPlane), 'Marker', SurfaceMarkerType(iBody), 'MarkerSize', GraphsMarkerSize)
        hold on
        grid on
        xlabel('Span [m]', 'FontSize', AxisLabelSize)
        ylabel('C_{L, 2D} [-]', 'FontSize', AxisLabelSize)
        title(sprintf('C_{L} Distribution, \\alpha = %.1f', Alpha(iAlpha)), '', 'FontSize', TitleSize)

    end 

    ax2 = subplot(1, 2, 2);
    for iBody = 1:config.NBodies % CD

        plot(CPoint{iBody,1}(:,2),CD_2Di{iBody,1},'DisplayName', sprintf('%s %s', PlaneName(iPlane), WingOrTail(iBody)), 'Color', PlaneLineColor(iPlane), 'LineWidth', GraphsLineWidth, 'LineStyle', SurfaceLineType(iBody), 'MarkerFaceColor', PlaneMarkerColor(iPlane), 'Marker', SurfaceMarkerType(iBody), 'MarkerSize', GraphsMarkerSize)
        hold on
        grid on
        xlabel('Span [m]', 'FontSize', AxisLabelSize)
        ylabel('C_{D, 2D} [-]', 'FontSize', AxisLabelSize)
        title(sprintf('C_{D} Distribution, \\alpha = %.1f', Alpha(iAlpha)), '', 'FontSize', TitleSize);

    end


end 

%% show all legends
for iSector = 1:FigCounter-1
    figure(iSector)
    legend('show')
    legend('FontSize', LegendSize)
end

legend(ax1, 'show', 'FontSize', LegendSize)
legend(ax2, 'show', 'FontSize', LegendSize)
    
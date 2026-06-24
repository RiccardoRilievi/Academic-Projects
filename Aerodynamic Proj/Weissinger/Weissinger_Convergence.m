close all
clear all
clc

Span = 80;
Alpha = 1;
Chord = 1;
U_inf = 1;
Rho = 1.225;

Discr_Prandtl = [1, 2, 3, 4, 5, 6, 7, 10, 20, 100];
Discr_Weissinger = [1, 5, 10, 15, 20, 25, 30, 40, 50, 60];

%% Prandtl 

for idiscr_p = 1: length(Discr_Prandtl)

    [L_Prandtl(idiscr_p), D_Prandtl(idiscr_p), CL_Prandtl(idiscr_p), CD_Prandtl(idiscr_p)] =  PrandtlFun(Span, Alpha, Chord, Discr_Prandtl(idiscr_p), U_inf, Rho);

end 

for idiscr_w = 1: length(Discr_Weissinger)
    %% Test Case 1
    
    U_Inf_Mag = U_inf;
    beta = 0;
    Alpha_U_inf = Alpha;
    U_Inf = [cosd(Alpha_U_inf)*cosd(beta) sind(beta) sind(Alpha_U_inf)] ./ norm([cosd(Alpha_U_inf) * cosd(beta) sind(beta) sind(Alpha_U_inf)]) .* U_Inf_Mag;
    rho = Rho;
    
    config.NBodies = 1;
    
    config.RootChord = [Chord];
    config.DihedralAngle = [0,0,0]; % [°]
    config.SweepAngle = [0,0,0];    % [°]
    config.TaperRatio = [1,1,1]; % TipChord/RootChord  
    config.Span = [Span];
    config.LEPosition_X = [0,10,0];
    config.LEPosition_Y = [0,0,0];
    config.LEPosition_Z = [0,0,0];
    
    config.RotationAngle_X = [0,0,0];   
    config.RotationAngle_Y = [0,0,0];   
    config.RotationAngle_Z = [0,0,0];   
    
    % Discretization options
    config.SemiSpanwiseDiscr = [Discr_Weissinger(idiscr_w),10,20];
    config.ChordwiseDiscr = [Discr_Weissinger(idiscr_w),10,20];
    
    
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
    
    ControlPoints = cell(config.NBodies, 1);
    InducedPoints = cell(config.NBodies, 1);
    Normals = cell(config.NBodies, 1);
    InfiniteVortices = cell(config.NBodies, 1);
    Vortices = cell(config.NBodies, 1);
    internalMesh = cell(config.NBodies, 1);
    WingExtremes = cell(config.NBodies, 1);
    
    % Create L, D
    
    L_2Di = cell(config.NBodies, 1);
    D_2Di = cell(config.NBodies, 1);
    CL_2Di = cell(config.NBodies, 1);
    CD_2Di = cell(config.NBodies, 1);
    U_Section = cell(config.NBodies, 1);
    L_3D = cell(config.NBodies, 1);
    CL_3D = cell(config.NBodies, 1);
    D_3D = cell(config.NBodies, 1);
    CD_3D = cell(config.NBodies, 1);
    
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
    
     %% Compute control point for each sector 
    CPoint = cell(config.NBodies, 1);
    
    for iBody = 1:config.NBodies
        CPoint{iBody, 1} = zeros(config.SemiSpanwiseDiscr(iBody), 3);
        for i = 1:config.SemiSpanwiseDiscr(iBody)
            CPoint_LE = (internalMesh{iBody,1}{1,i}.LERoot + internalMesh{iBody,1}{1,i}.LEtip)/2;
            CPoint_TE = (internalMesh{iBody,1}{end,i}.TERoot + internalMesh{iBody,1}{end,i}.TEtip)/2;
        
            versor = (CPoint_LE - CPoint_TE) / norm(CPoint_LE - CPoint_TE);
        
            CPoint{iBody, 1}(i,:) = CPoint_TE + (3/4 * norm(CPoint_LE - CPoint_TE) ) * versor;
            
        end 
        CPoint{iBody, 1} = [CPoint{iBody, 1}; flip([CPoint{iBody, 1}(:,1), -CPoint{iBody, 1}(:,2), CPoint{iBody, 1}(:,3)])];
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
    
    U_S = [];
    U_Section_s = zeros(1,3);
    
    
    
    for iBody = 1:config.NBodies
            % Cycle on sectors
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
                            [U_Section_s] = U_Section_s + Gamma{jCorpo}(ChordPanel_j, SpanPanel_j)*vortexInfluence(CPoint{iBody,1}(Sectors_i,:), Extreme_1, Extreme_2);
    
    
    
                            % Compute the influence induced by second
                            % semi-infinite vortex
                            Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.onWing;
                            Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.toInfty;
                            [U_Section_s] = U_Section_s + Gamma{jCorpo}(ChordPanel_j, SpanPanel_j)*vortexInfluence(CPoint{iBody,1}(Sectors_i,:), Extreme_1, Extreme_2);
    
    
    
                        end
                    end
    
                end
    
               U_S = [U_S; U_Section_s];
               U_Section_s = zeros(1,3);
    
            end
    
             U_Section{iBody,1} = U_S;
             clear U_S
             U_S = [];
    end
    
    %% Compute the 2D and 3D Lift
    
    L_2Di_s = [];
    
    
    for iBody = 1:config.NBodies
        for i = 1:config.SemiSpanwiseDiscr(iBody)*2
            L2Di = sum(rho * U_Inf_Mag * Gamma{iBody}(:,i)*cosd(config.DihedralAngle(iBody)));
            L_2Di_s=[L_2Di_s, L2Di];
            L_2Di{iBody}= L_2Di_s;
        end 
        clear L_2Di_s
        L_2Di_s = [];
    
        CL_2Di{iBody} = L_2Di{iBody}./( 0.5*rho*U_Inf_Mag^2 *config.Surface(iBody));
    end
    
    
    Delta_bi = cell(config.NBodies, 1);
    for iBody = 1:config.NBodies
        Delta_bi{iBody, 1} = (config.SemiSpan(iBody)*2 / (config.SemiSpanwiseDiscr(iBody)*2)) * ones(config.SemiSpanwiseDiscr(iBody)*2, 1);
    
        L_3D{iBody, 1} = sum(Delta_bi{iBody, 1}.*L_2Di{iBody, 1}');
        CL_3D{iBody, 1} = L_3D{iBody, 1}./( 0.5*rho*U_Inf_Mag^2 *config.Surface(iBody));
    
    end
         
    
    
    
    %% Compute 2D and 3D induced drag
    
    alpha_i = cell(config.NBodies, 1);
    
    for iBody = 1:config.NBodies
        for i=1:config.SemiSpanwiseDiscr(iBody)*2
        
                alphai = atan2( (dot(U_Section{iBody,1}(i,:),Normals{iBody,1}{1,i}.Coords)), U_Inf_Mag);
                alpha_i{iBody,1} = [alpha_i{iBody,1}; alphai];
         
        end 
    
        D_2Di{iBody,1} = -L_2Di{iBody,1}' .* sin((alpha_i{iBody,1}));
        CD_2Di{iBody,1} = D_2Di{iBody,1}./( 0.5*rho*U_Inf_Mag^2 *config.Surface(iBody));
    
        D_3D{iBody,1} = sum(D_2Di{iBody,1}.*Delta_bi{iBody,1});
        CD_3D{iBody,1} = D_3D{iBody,1}./( 0.5*rho*U_Inf_Mag^2 *config.Surface(iBody));
    
    end 

    for iBody = 1:config.NBodies
            CL_Weissinger(idiscr_w,iBody) = CL_3D{iBody, 1};
            CD_Weissinger(idiscr_w,iBody) = CD_3D{iBody, 1};
            L_Weissinger(idiscr_w,iBody) = L_3D{iBody, 1};
            D_Weissinger(idiscr_w,iBody) = D_3D{iBody, 1};
    end 

end 

%% Graph 
close all
subplot(1,2,1)
plot(Discr_Weissinger, L_Weissinger,'-o','LineWidth',2)
hold on
grid on
plot(Discr_Prandtl, L_Prandtl,'-o','LineWidth',2)
xlabel('Discretization','FontSize', 20)
ylabel('Lift [N]','FontSize', 20)
legend('Weissinger', 'Prandtl', 'FontSize', 15)
title('Convergence','FontSize', 20)

subplot(1,2,2)
plot(Discr_Weissinger, D_Weissinger,'-o','LineWidth',2)
hold on
grid on
plot(Discr_Prandtl, D_Prandtl,'-o','LineWidth',2)
xlabel('Discretization','FontSize', 20)
ylabel('Drag [N]','FontSize', 20)
legend('Weissinger', 'Prandtl', 'FontSize', 15)
title('Convergence','FontSize', 20)

T_Weissinger = table( ...
            Discr_Weissinger(:), ...
            CL_Weissinger(:,:), ...
            CD_Weissinger(:,:), ...
            L_Weissinger(:,:), ...
            D_Weissinger(:,:), ...
            'VariableNames', {'Discretization','CL','CD','L','D'} );

disp(T_Weissinger)

T_Prandtl = table( ...
    Discr_Prandtl(:), ...
    CL_Prandtl(:), ...
    CD_Prandtl(:), ...
    L_Prandtl(:), ...
    D_Prandtl(:), ...
    'VariableNames', {'Discretization','CL','CD','L','D'} );

disp(T_Prandtl)

T_delta = table( ...
            Discr_Weissinger(:),...
            abs(CL_Weissinger(:)-CL_Prandtl(:)),...
            abs(CD_Weissinger(:)-CD_Prandtl(:)),...
            abs(L_Weissinger(:)-L_Prandtl(:)),...
            abs(D_Weissinger(:)-D_Prandtl(:)),...
            'VariableNames', {'Discretization','\DeltaCL','\DeltaCD','\DeltaL','\DeltaD'} );

disp(T_delta)

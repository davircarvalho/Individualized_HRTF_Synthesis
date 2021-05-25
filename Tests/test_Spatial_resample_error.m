clear all; clc; 
% DAVI ROCHA CARVALHO; ENG. ACUSTICA - UFSM; Abril/2021
%% PATHs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath([pwd, '\..\Functions']);
addpath([pwd, '\..\DADOS_TREINAMENTO']);
local = [pwd, '\..\Datasets\'];

% CIPIC
pathcipic = dir([local 'CIPIC\*.sofa']);
[~,idx_cipic] = natsortfiles({pathcipic.name});
pathcipic = pathcipic(idx_cipic, :);

% ARI
pathari = dir([local 'ARI\hrtf b_nh*.sofa']);
[~,idx_ari] = natsortfiles({pathari.name});
pathari = pathari(idx_ari, :);

% ITA
pathita = dir([local 'AACHEN\*.sofa']);
[~,idx_ita] = natsortfiles({pathita.name});
pathita =  pathita(idx_ita, :);

% 3D3A
path3d3a = dir([local '3D3A\Public-Data\Subject*\Subject*_HRIRs.sofa']);
[~,idx_3d3a] = natsortfiles({path3d3a.name});
path3d3a =  path3d3a(idx_3d3a, :);

% % TU Berlim 
pathtub_meas = dir([local 'HUTUBS\pp*_HRIRs_measured.sofa']);
[~,idx_tubmeas] = natsortfiles({pathtub_meas.name});
pathtub_meas = pathtub_meas(idx_tubmeas, :);

% TU Berlim 
pathtub_sim = dir([local 'HUTUBS\pp*_HRIRs_simulated.sofa']);
[~,idx_tubsim] = natsortfiles({pathtub_sim.name});
pathtub_sim = pathtub_sim(idx_tubsim, :);

% Viking
pathtub_vik = dir([local 'VIKING\*.sofa']);
[~,idx_vik] = natsortfiles({pathtub_vik.name});
pathtub_vik = pathtub_vik(idx_vik, :);



k = 10;
CIPIC = SOFAload([pathcipic(k).folder '\' pathcipic(k).name]);
ARI = SOFAload([pathari(k).folder, '\',pathari(k).name]);
ITA = SOFAload([pathita(k).folder, '\',pathita(k).name]);
TUBmeas = SOFAload([pathtub_meas(k).folder '\' pathtub_meas(k).name]);   
TUBsim = SOFAload([pathtub_sim(k).folder '\' pathtub_sim(k).name]);                

D3A = SOFAload([path3d3a(k).folder '\' path3d3a(k).name]);       
VIK = SOFAload([pathtub_vik(k).folder '\' pathtub_vik(k).name]);       


%%% Transição de coordenadas cartesianas para esfericas
for l = 1:length(ITA.SourcePosition)
    x = ITA.SourcePosition(l, 1);  
    y = ITA.SourcePosition(l, 2); 
    z = ITA.SourcePosition(l, 3);
    % new coordinates
    [az,elev,r] = cart2sph(x,y,z);
    azi=rad2deg(az); elev=rad2deg(elev);
    [azi,ele] = nav2sph(azi,elev);
    azi(azi == 360) = 0;
    % update coordinates
    ITA.SourcePosition(l, :) = [azi, ele, round(r)];
    % more metadata
    ITA.SourcePosition_Type = 'spherical';
    ITA.SourcePosition_Units = 'degree, degree, meter';              
end       



%% Plot erro
out_pos = D3A.SourcePosition;
fitgrid(CIPIC, out_pos, 'CIPIC')
fitgrid(ARI, out_pos, 'ARI')
fitgrid(ITA, out_pos, 'AACHEN')
fitgrid(TUBsim, out_pos, 'HUTUBS')
% fitgrid(VIK, out_pos, ' VIKING')


function fitgrid(Obj, out_pos, titl)
    [~, erro] = sofaFit2Grid(Obj, out_pos, 'adapt');
    h = figure();
    scatter(out_pos(:,1), out_pos(:,2), 25, erro, 'filled')
    % properties 
    c = colorbar;
    c.Label.String = 'RMS error (°)';
    caxis([0 10]); 
    title(titl)
    axis tight
    xlabel('Azimuth (°)')
    ylabel('Elevation (°)')
    set(gca, 'FontSize', 13)
    filename=['posi_error_' titl '.pdf'];
    exportgraphics(h,filename,'ContentType','vector')
end

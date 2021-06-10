% clear all; clc; tic
% DAVI ROCHA CARVALHO; ENG. ACUSTICA - UFSM; DEZEMBRO/2020
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

% RIEC
pathriec = dir([local 'RIEC\*.sofa']);
[~,idx_riec] = natsortfiles({pathriec.name});
pathriec =  pathriec(idx_riec, :);

% TU Berlim 
pathtub_meas = dir([local 'HUTUBS\pp*_HRIRs_measured.sofa']);
[~,idx_tubmeas] = natsortfiles({pathtub_meas.name});
pathtub_meas = pathtub_meas(idx_tubmeas, :);


pathtub_sim = dir([local 'HUTUBS\pp*_HRIRs_simulated.sofa']);
[~,idx_tubsim] = natsortfiles({pathtub_sim.name});
pathtub_sim = pathtub_sim(idx_tubsim, :);

% VIKING
pathvik = dir([local 'VIKING\*.sofa']);
[~,idx_vik] = natsortfiles({pathvik.name}); % garantir que estão em ordem
pathvik = pathvik(idx_vik, :);


%% Output positions
res = 1;
max_el = 90;
min_el = -90;
elevation=[min_el:res:max_el max_el-1:-res:min_el zeros(1,length(1:res:360-res))]';
azimuth=[zeros(length(min_el:res:max_el),1); 180*ones(length(max_el-1:-res:min_el),1); (1:res:360-res)'];
     
out_pos = [azimuth, elevation];

%%
k = 5;
CIPIC = SOFAload([pathcipic(k).folder '\' pathcipic(k).name]);
ARI = SOFAload([pathari(k).folder, '\',pathari(k).name], 'nochecks');
ITA = SOFAload([pathita(k).folder, '\',pathita(10).name], 'nochecks');
D3A = SOFAload([path3d3a(k).folder '\' path3d3a(k).name], 'nochecks');       
TUBmeas = SOFAload([pathtub_meas(k).folder '\' pathtub_meas(k).name], 'nochecks');                
TUBsim = SOFAload([pathtub_sim(k).folder '\' pathtub_sim(k).name], 'nochecks');
VIK = SOFAload([pathvik(k).folder '\' pathvik(k).name], 'nochecks');

%%
% % % [ITA] Transição de coordenadas cartesianas para esfericas
for l = 1:length(ITA.SourcePosition)
    x = ITA.SourcePosition(l, 1);  
    y = ITA.SourcePosition(l, 2); 
    z = ITA.SourcePosition(l, 3);
    % new coordinates
    [az,elev,r] = cart2sph(x,y,z);
    azi=rad2deg(az); elev=rad2deg(elev);
    [azi,ele]   = nav2sph(azi,elev);
    % update coordinates
    ITA.SourcePosition(l, 1) = azi;
    ITA.SourcePosition(l, 2) = ele; 
    ITA.SourcePosition(l, 3) = round(r);
    % more metadata
    ITA.SourcePosition_Type = 'spherical';
    ITA.SourcePosition_Units = 'degree, degree, meter';              
end       


%%
close all; clc
% CIPIC_interp = process2unite(CIPIC, out_pos, 'cipic');
% ARI_interp = process2unite(ARI, out_pos,  'ari');
% ITA_interp = process2unite(ITA, out_pos,  'ita');
D3A_interp = process2unite(D3A, out_pos,  'd3a');
% TUBmeas_interp = process2unite(TUBmeas, out_pos,  'tubmeas');
% TUBsim_interp = process2unite(TUBsim, out_pos, 'tubsim');
% VIK_interp = process2unite(VIK, out_pos, 'VIKING');


function Obj_out = process2unite(Obj, out_pos, name)
    % Make same grid
    Obj_out = SOFA_barycentric_interp(Obj, out_pos);     

    
    %% PLOTA
    plane1 = 'MagMedian';
    plane2 = 'MagHorizontal';
    
    figure()
    SOFAplotHRTF(Obj,plane1); title(['Reference - ' plane1, ' ' name]);
    axis tight
    xlim([0 2e4])
        
    figure()
    SOFAplotHRTF(Obj,plane2); title(['Reference - ' plane2, ' ' name]);
    axis tight
    xlim([0 2e4])
     
    figure()
    SOFAplotHRTF(Obj_out,plane1); title(['Interpolated - ' plane1, ' ' name]);
    axis tight
    xlim([0 2e4])
    
    figure()
    SOFAplotHRTF(Obj_out,plane2); title(['Interpolated - ' plane2, ' ' name]);
    axis tight
    xlim([0 2e4])
    pause(0)
end

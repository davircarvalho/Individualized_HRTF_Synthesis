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


%% Options
% Defina quais datasets usar: {'cipic', 'ari', 'ita', '3d3a', 'riec', 'tub_meas', 'tub_sim'}, o
Datasets = {'cipic', 'ari', 'ita', '3d3a', 'tub_meas', 'tub_sim'};
no_samples = 200; % Tamanho do vetor na saída (pós fft)
fs   = 44100;     % Taxa de amostragem 
fmin = 250;       % Frequencia min de corte para RI 
fmax = 18000;     % Frequencia max de corte para IR  
% Grid objetivo selecionado a partir do grid com menor número de posições
% out_pos = select_best_grid(Datasets);
freq = linspace(0, fs-fs/no_samples, no_samples);


res = .1;
ele=[-90:res:90 89:-res:-90 zeros(1,length(1:res:360-res))]';
azi=[zeros(length(-90:res:90),1); 180*ones(length(89:-res:-90),1); (1:res:360-res)'];

out_pos = [azi, ele, ones(length(azi), 2)];


%%
k = 10;
% CIPIC = SOFAload([pathcipic(k).folder '\' pathcipic(k).name]);
% ARI = SOFAload([pathari(k).folder, '\',pathari(k).name], 'nochecks');
% ITA = SOFAload([pathita(k).folder, '\',pathita(10).name], 'nochecks');
% D3A = SOFAload([path3d3a(k).folder '\' path3d3a(k).name], 'nochecks');       
TUBmeas = SOFAload([pathtub_meas(k).folder '\' pathtub_meas(k).name], 'nochecks');                
TUBsim = SOFAload([pathtub_sim(k).folder '\' pathtub_sim(k).name], 'nochecks');

%%
% [ITA] Transição de coordenadas cartesianas para esfericas
% for l = 1:length(ITA.SourcePosition)
%     x = ITA.SourcePosition(l, 1);  
%     y = ITA.SourcePosition(l, 2); 
%     z = ITA.SourcePosition(l, 3);
%     % new coordinates
%     [az,elev,r] = cart2sph(x,y,z);
%     azi=rad2deg(az); elev=rad2deg(elev);
%     [azi,ele]   = nav2sph(azi,elev);
%     % update coordinates
%     ITA.SourcePosition(l, 1) = azi;
%     ITA.SourcePosition(l, 2) = ele; 
%     ITA.SourcePosition(l, 3) = round(r);
%     % more metadata
%     ITA.SourcePosition_Type = 'spherical';
%     ITA.SourcePosition_Units = 'degree, degree, meter';              
% end       


%%
close all
% CIPIC_interp = process2unite(CIPIC, out_pos, fs, 'cipic');
% ARI_interp = process2unite(ARI, out_pos, fs, 'ari');
% ITA_interp = process2unite(ITA, out_pos, fs, 'ita');
% D3A_interp = process2unite(D3A, out_pos, fs, 'd3a');
TUBmeas_interp = process2unite(TUBmeas, out_pos, fs, 'tubmeas');
% TUBsim_interp = process2unite(TUBsim, out_pos, fs, 'tubsim');










function Obj_out = process2unite(Obj, out_pos, fs, name)
    % Make same grid
    Obj_out = sofaFit2Grid(Obj, out_pos, 'spherical_harmonics', 'Fs', fs);     
    % Normalize L/R balance and IR levels
%     Obj = sofaNormalize(Obj);
%     % filter
%     Obj = sofaIRfilter(Obj, fmin, fmax);
%     % HRTF -> DTF
%     [Obj, ~] = SOFAhrtf2dtf(Obj);    



    %% PLOTA
    plane1 = 'MagHorizontal';
    plane2 = 'MagSagittal';
    
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
clear all; clc; tic
% DAVI ROCHA CARVALHO; ENG. ACUSTICA - UFSM; DEZEMBRO/2020
%% PATHs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addpath([pwd, '\..\Functions']);
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
% res = 1;
% ele=[-90:res:90 89:-res:-90 zeros(1,length(1:res:360-res))]';
% azi=[zeros(length(-90:res:90),1); 180*ones(length(89:-res:-90),1); (1:res:360-res)'];
% r = ones(length(azi),1);
% 
% out_pos = unique([azi, ele, r], 'rows');
out_pos = equiangular_coordinates(4);


%%
k = 1;
CIPIC = SOFAload([pathcipic(k).folder '\' pathcipic(k).name]);
ARI = SOFAload([pathari(k).folder, '\',pathari(k).name], 'nochecks');
ITA = ITA2spheric(SOFAload([pathita(k).folder, '\',pathita(k).name], 'nochecks'));
D3A = SOFAload([path3d3a(k).folder '\' path3d3a(k).name], 'nochecks');       
TUBmeas = SOFAload([pathtub_meas(k).folder '\' pathtub_meas(k).name], 'nochecks');                
TUBsim = SOFAload([pathtub_sim(k).folder '\' pathtub_sim(k).name], 'nochecks');
VIK = SOFAload([pathvik(k).folder '\' pathvik(k).name], 'nochecks');

% TUBsim.SourcePosition = inpt_pos;
% TUBsim.Data.IR = inpt_IR;
%%
close all
CIPIC_interp = process2unite(CIPIC, out_pos,  'cipic');
% ARI_interp = process2unite(ARI, out_pos,  'ari');
% ITA_interp = process2unite(ITA, out_pos, 'ita'); 
% D3A_interp = process2unite(D3A, out_pos, 'd3a');
% TUBmeas_interp = process2unite(TUBmeas, out_pos, 'tubmeas');
% TUBsim_interp = process2unite(TUBsim, out_pos, 'tubsim');
% VIK_interp = process2unite(VIK, out_pos, 'VIKING');


% SOFAplotGeometry(ITA_interp)


function Obj_out = process2unite(Obj, out_pos, name)
    % Make same grid
    Obj_out = sofaSHinterpolate(Obj, out_pos, 'ITA_api');
%     Obj_out = sofaFit2Grid(Obj, out_pos, 'spherical_harmonics');     

% Check sd 
    fmin = 500; fmax = 18000;
    SD = mean(sofaSpecDist(Obj_out, Obj, fmin,fmax), 'all')
    %% PLOTA
    plane1 = 'MagHorizontal';
    plane2 = 'MagSagittal';
    
%     figure()
%     SOFAplotHRTF(Obj,plane1); title(['Reference - ' plane1, ' ' name]);
%     xlim([0 2e4])
%     ylim([-180, 180])
%     
%     figure()
%     SOFAplotHRTF(Obj,plane2); title(['Reference - ' plane2, ' ' name]);
%     xlim([0 2e4])
%     ylim([-80, 250])
%      
    figure()
    SOFAplotHRTF(Obj_out,plane1); title(['Interpolated - ' plane1, ' ' name]);
    xlim([0 2e4])
    ylim([-180, 180])
    
    figure()
    SOFAplotHRTF(Obj_out,plane2); title(['Interpolated - ' plane2, ' ' name]);
    xlim([0 2e4])
    ylim([-80, 250])
end
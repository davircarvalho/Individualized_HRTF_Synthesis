clear all; clc; 
% DAVI ROCHA CARVALHO; ENG. ACUSTICA - UFSM; NOVEMBRO/2020
% Calculo do ITD para HRIRs Fabian 
%% Load HRIRs
path = [pwd '\..\Datasets\Generic HRTFs\'];
Fabian = SOFAload([path 'FABIAN_HRIR_measured_HATO_0.sofa']);

%% Select desired positions
addpath(genpath([pwd, '\..\DADOS_TREINAMENTO'])); 
addpath(genpath([pwd, '\..\Functions'])); 
load('DADOS_TREINAMENTO\target_pca_CIPIC_ARI_ITA_3D3A.mat');
Fabian = sofaFit2Grid(Fabian, out_pos, 'spherical_harmonics');

%% ITD estimate
ref_itd = sofaGetITD(Fabian, 'time', 'thr', 20);

figure()
plot(ref_itd)


%% More data
ref_pos = Fabian.SourcePosition;

% antropometria
% fabian_hwidth = 15.5435; %[cm]
% fabian_hdepth = 19.9160; %[cm]
ref_width = 15.497; %[cm]
ref_depth = 19.750; %[cm]
ref_height = 24;     %[cm]


%% Save
save('..\Functions\fabian_itd.mat', 'ref_itd', 'ref_pos', 'ref_width', ...
                       'ref_depth', 'ref_height')
disp('Dados Salvos!')

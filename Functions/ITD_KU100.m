clear all; clc; 
% DAVI ROCHA CARVALHO; ENG. ACUSTICA - UFSM; MARÇO/2020
% Calculo do ITD para HRIRs Fabian 
%% Load HRIRs
Obj = SOFAload('ClubFritz5.sofa'); %medição do ITA

%% Desired ITD positions
addpath(genpath([pwd, '\..\DADOS_TREINAMENTO'])); 
load('DADOS_TREINAMENTO\target_pca_CIPIC_ARI_ITA_3D3A.mat');
Obj = sofaFit2Grid(Obj, out_pos, 'spherical_harmonics');


%% ITD estimate
ref_itd = sofaGetITD(Obj, 'time'); % use: ('correlation', true)
figure()
plot(ref_itd)

%% More data
ref_pos = Obj.SourcePosition;

% antropometria
ref_width  = 15.5; %[cm]
ref_depth  = 20;   %[cm]
ref_height = 25;   %[cm]

%% Save
save('Functions\KU100_itd.mat', 'ref_itd', 'ref_pos', 'ref_width', ...
                       'ref_depth', 'ref_height')
disp('Dados Salvos!')

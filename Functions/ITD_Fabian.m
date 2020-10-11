clear all; clc; 
% DAVI ROCHA CARVALHO; ENG. ACUSTICA - UFSM; MARÇO/2020
% Calculo do ITD para HRIRs Fabian 
%% Load HRIRs
path = 'B:\Documentos\#3 - TCC\EAC-TCC-Davi\HRTF-datasets\Cabecas\';
Fabian = SOFAload([path 'FABIAN_HRIR_measured_HATO_0.sofa']);

%% ITD estimate
ref_itd = sofaGetITD(Fabian, 'time');

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
save('Functions\fabian_itd.mat', 'ref_itd', 'ref_pos', 'ref_width', ...
                       'ref_depth', 'ref_height')
disp('Dados Salvos!')

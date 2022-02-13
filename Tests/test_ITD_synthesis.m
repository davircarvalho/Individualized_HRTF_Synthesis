clear all; 


%% LOAD TEST DATA 
% Load Antropometria 
load('..\DADOS_TREINAMENTO\input_CIPIC.mat');

%%% HRIRs ORIGINAIS
local = [fileparts(which('test_ITD_synthesis')), '\..\Datasets\'];
pathcipic = dir([local 'HUTUBS\pp*_HRIRs_simulated.sofa']);
subj = 3;
Obj = SOFAload([pathcipic(subj).folder '\' pathcipic(subj).name], 'nochecks');

Obj.Data.IR = Obj.Data.IR./max(abs(Obj.Data.IR(:)));
%% Properties
pos = Obj.SourcePosition;
fs = Obj.Data.SamplingRate;
width = anthro(1,subj,1);
depth = anthro(2,subj,1);


%% calcular ITD
method = 'spheric';
itd_sphe = itd_synthesis(width, depth, pos, fs, method, 'time');

method = 'adapt';
itd_adpt = itd_synthesis(width, depth, pos, fs, method, 'time');

itd_real = SOFAgetITD(Obj, 'time', 'thr', 20, 'upsample', 32);

sum(abs(itd_real - itd_sphe))
sum(abs(itd_real - itd_adpt))

%% PLOT
figure
plot(itd_real, 'linewidth', 1.5); hold on
plot(itd_sphe);
plot(itd_adpt); hold off
legend('Measured', 'Spheric', 'Adapt', 'location', 'best')
ylabel('Time (s)')
xlabel('Position index')
axis tight

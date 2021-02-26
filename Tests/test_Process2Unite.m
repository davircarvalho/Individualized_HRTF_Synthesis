clear all; clc
addpath(genpath([pwd '\..\Functions']))
%% Reference
local = 'D:\Documentos\1 - Work\Individualized_HRTF_Synthesis\Datasets\HUTUBS\';
Obji = SOFAload([local 'pp1_HRIRs_measured.sofa']);
[Obj, ~] = SOFAhrtf2dtf(Obji); 
idx = dsearchn(Obj.SourcePosition(:,[1,2]), [0,0]);
IR = squeeze(Obj.Data.IR(idx,1,:));
N = length(IR);
fs = Obj.Data.SamplingRate;
f = linspace(0,fs-fs/N, N);
    
%% Parameters
out_pos = Obj.SourcePosition;
fmin = 250;
fmax = 18000;

%% Process2Unite
Obj_u = process2unite(Obji, out_pos, fs, fmin, fmax);
idx_u = dsearchn(Obj_u.SourcePosition(:,[1,2]), [0,0]);
IR_u = squeeze(Obj_u.Data.IR(idx_u,1,:));
N_u = length(IR_u);
fs_u = Obj_u.Data.SamplingRate;
f_u = linspace(0,fs_u-fs_u/N_u, N_u);

%% plot
% close all;
figure()
plot(f, db(abs(fft(IR(:,1)))));hold on 
plot(f_u, db(abs(fft(IR_u(:,1)))));hold off
xlim([0, 2.4e4])


%% LOCAL FUCTIONS 
function Obj = process2unite(Obj, out_pos, fs, fmin, fmax)
    % Make same grid
    Obj = sofaFit2Grid(Obj, out_pos, 'adapt', 'Fs', fs);     
    % Normalize L/R balance and IR levels
    Obj = sofaNormalize(Obj);
    % filter
    Obj = sofaIRfilter(Obj, fmin, fmax);
    % HRTF -> DTF
    [Obj, ~] = SOFAhrtf2dtf(Obj); 
end
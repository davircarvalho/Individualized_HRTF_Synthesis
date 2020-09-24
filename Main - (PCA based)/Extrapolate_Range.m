% Extrappolação de hrtf quanto ao raio 
% Davi R Carvalho MAIO/2020
clear all; clc; tic
%% LOAD  
Obj = SOFAload('individuo_140.sofa');
Obj = sofaResample(Obj, 48000);
fs = Obj.Data.SamplingRate;
NeqDataset = 44;

%% Converter para SUpDEq
HRIRs_sfd = supdeq_sofa2sfd(Obj, NeqDataset, [], 2); 
sparseHRTFdataset = supdeq_sofa2hrtf(Obj, NeqDataset, [], 2);

%% GRID  
% out_pos = sparseHRTFdataset.samplingGrid;
% fix coordinates to SUpDEq format
out_pos = Obj.SourcePosition;
out_pos(:,2) = 90-out_pos(:,2);


%% Equalize The sphere with the input HRTFs
%The eqDataset describes the sound pressure distribution on a sphere 
%Use N = 44 and defaults: earDistance = 0.165m, NFFT = 512, fs = 48000;
eqDataset = supdeq_getEqDataset(NeqDataset, [],...
                                (2*size(sparseHRTFdataset.HRTF_L, 2)-1),...
                                fs);
                            
eqHRTFdataset = supdeq_eq(sparseHRTFdataset,eqDataset,NeqDataset, out_pos);


%% deEQUALIZER
targetDistance = 0.3; %[m] nova distancia
deqDataset = supdeq_getEqDataset(NeqDataset,[],[],fs,1,targetDistance);

% HRIRs_sfd_ext = supdeq_rangeExt(HRIRs_sfd, targetDistance);
% deqDataset.Hl_nm = HRIRs_sfd_ext.Hl_nm;
% deqDataset.Hr_nm = HRIRs_sfd_ext.Hr_nm;

[denseHRTFdataset, denseHRIRdataset, denseHRTFdataset_sh] = supdeq_deq(eqHRTFdataset,...
                                                                       deqDataset,...
                                                                       NeqDataset,...
                                                                       out_pos,...
                                                                       [8,32],0.99);

%% SOFA output
Obj_out = supdeq_writeSOFAobj(denseHRIRdataset.HRIR_L,...
                              denseHRIRdataset.HRIR_R,...
                              out_pos,fs,[],...
                              targetDistance);

file_path = fullfile(['supdeq_5m.sofa']);
compression = 0;
SOFAsave(file_path, Obj_out, compression);
toc
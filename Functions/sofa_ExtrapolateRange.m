function Obj_out = sofa_ExtrapolateRange(Obj, targetDistance)
% Extrappolação de hrtf quanto ao raio 
% Davi R Carvalho MAIO/2020
%% LOAD  
% Obj = SOFAload('individuo_140.sofa');
% fs = Obj.Data.SamplingRate;
% out_pos(:,2) = 90-out_pos(:,2);
% out_pos(:,1) = 360-out_pos(:,1);
% out_pos(out_pos == 360, 1) = 0;
fs = Obj.SamplingRate;
%% Converter para SUpDEq
% HRIRs_sfd = supdeq_sofa2sfd(Obj,35,[],2); 
% sparseHRTFdataset = supdeq_sofa2hrtf(Obj,4, [],2);
sparseHRTFdataset = supdeq_sofa2hrtf(Obj,4, [], 2);
sparseSamplingGrid = sparseHRTFdataset.samplingGrid;
out_pos = sparseSamplingGrid;
%% Equalize The sphere with the input HRTFs
%The eqDataset describes the sound pressure distribution on a sphere 
%Use N = 44 and defaults: earDistance = 0.165m, NFFT = 512, fs = 48000;
NeqDataset = 44;
eqDataset = supdeq_getEqDataset(NeqDataset,[], (2*size(sparseHRTFdataset.HRTF_L, 2)-1), fs);
% sparseSamplingGrid = out_pos;
Nsparse = sparseHRTFdataset.Nmax;
eqHRTFdataset = supdeq_eq(sparseHRTFdataset,eqDataset,Nsparse, sparseSamplingGrid);


%%
% targetDistance = 1; %[m] nova distancia
Ndense = 44;
deqDataset = supdeq_getEqDataset(NeqDataset,[],[],[],1,targetDistance);
denseSamplingGrid = out_pos;

[~, denseHRIRdataset, ~] = supdeq_deq(eqHRTFdataset,...
                                    deqDataset, Ndense, denseSamplingGrid);


Obj_out = supdeq_writeSOFAobj(denseHRIRdataset.HRIR_L,...
                                denseHRIRdataset.HRIR_R,...
                                denseSamplingGrid,[],[],...
                                targetDistance);

% file_path = fullfile(['supdeq_5m.sofa']);
% compression = 0;
% SOFAsave(file_path, Obj_out, compression);
% toc


end
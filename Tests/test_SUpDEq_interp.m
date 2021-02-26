%Load sparse HRIR dataset in SOFA format
clear all; clc
local = [pwd '\..\Datasets\'];    
path = dir([local 'HUTUBS\pp*_HRIRs_simulated.sofa']);
Obj = SOFAload([path(1).folder, '\',path(1).name], 'nochecks');
sparseHRIRdataset_SOFA = Obj;

% sparseHRIRdataset_SOFA = SOFAload('sparseHRIRdataset_L38.sofa');
%Transform to sparseHRTFdataset struct with pre-defined samplingGrid 
%(Lebedev grid with 38 nodes here), Nmax = 4, and FFToversize = 4.
Nmax = floor(sqrt(Obj.API.M/4));
FFToversize = 4;
sparseHRTFdataset = supdeq_sofa2hrtf(sparseHRIRdataset_SOFA,Nmax,[],FFToversize);

%% (3) - Get equalization dataset (SH-coefficients)
%The eqDataset describes the sound pressure distribution on a sphere 
%Use defaults: N = 35, earDistance = 0.165m, NFFT = 512, fs = 48000;
NFFT=(length(sparseHRTFdataset.HRTF_L(1,:))-1)*2;
eqDataset = supdeq_getEqDataset(35, [], NFFT, Obj.Data.SamplingRate);


%% (4) - Perform equalization
%Here, the sparse HRTF dataset is equalized with the eqDataset. The
%equalized HRTF are transformed to the SH-domain again with the maximal 
%order N which is possible with the sparse sampling grid.
%N and the sparse sampling grid are part of the sparseHRTFdataset struct
sparseSamplingGrid = sparseHRTFdataset.samplingGrid;
Nsparse = sparseHRTFdataset.Nmax;

eqHRTFdataset = supdeq_eq(sparseHRTFdataset,eqDataset,Nsparse,sparseSamplingGrid,...
                            'tikhEps', 10^-4);


% (5) - Perform de-equalization 
%Here, the sparse equalized HRTF dataset is de-equalized with the
%deqDataset. This is done on a dense spatial sampling grid. The results is a
%dense HRTF/HRIR dataset. In this case, deqDataset is the same as the
%eqDataset...

%First, define dense spatial sampling grid. Here, we use the lebedev grid
%with 2702 points again (same as the reference HRIR dataset).
%The highest stable grid order here is N = 44. However, we use N = 35 for the
%spherical Fourier transform of the de-equalized HRTFs.
res = 1;
ele=[-90:res:90 89:-res:-90 zeros(1,length(1:res:355))]';
azi=[zeros(length(-90:res:90),1); 180*ones(length(89:-res:-90),1); (1:res:355)'];

des_pos = [azi, ele]; % sofa style 
des_pos(:,2) = 90-des_pos(:,2); % supdeq style


denseSamplingGrid = des_pos;
% denseSamplingGrid = supdeq_lebedev(2702);
Ndense = 35;

%Perform de-equalization. Apply head and tail window (8 and 32 samples
%respectively) to de-equalized HRIRs/HRTFs.
[~, denseHRIRdataset, ~] = ...
    supdeq_deq(eqHRTFdataset, eqDataset, Ndense, denseSamplingGrid);


% (6) - Optional: Save as SOFA object
%Use defaults: fs = 48000, earDistance = 0.165m, sourceDistance = 3.0m
denseHRIRdataset_SOFA = supdeq_writeSOFAobj(denseHRIRdataset.HRIR_L,...
                                            denseHRIRdataset.HRIR_R,...
                                            denseSamplingGrid);
                                        

%%
plane = 'MagSagittal';
close all
figure()
SOFAplotHRTF(sparseHRIRdataset_SOFA,plane); title(['Reference - ' plane]);
axis tight
figure()
SOFAplotHRTF(denseHRIRdataset_SOFA,plane); title(['Interpolated - ' plane]);
axis tight
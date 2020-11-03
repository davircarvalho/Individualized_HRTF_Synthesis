% SOFA HRTF INTERPOLATATION USING SPHERICAL HARMONICS
% function  Obj = sofaSHinterpolate(Obj, pos)
% ele = pos(:,2);
% azi = pos(:,1);
close all; clear all; clc

%% Import
% IR = Obj;
IR=SOFAload('D:\Documentos\1 - Work\Individualized_HRTF_Synthesis\Datasets\ARI\hrtf b_nh16.sofa');
fs=IR.Data.SamplingRate;
IR.GLOBAL_APIVersion=SOFAgetVersion;
%% Convert to TF
TF=SOFAgetConventions('SimpleFreeFieldHRTF');
TF.ListenerPosition=IR.ListenerPosition;
TF.ListenerPosition_Type=IR.ListenerPosition_Type;
TF.ListenerPosition_Units=IR.ListenerPosition_Units;
TF.ListenerView=IR.ListenerView;
TF.ListenerView_Type=IR.ListenerView_Type;
TF.ListenerView_Units=IR.ListenerView_Units;
TF.ListenerUp=IR.ListenerUp;
TF.SourcePosition=IR.SourcePosition;
TF.SourcePosition_Type=IR.SourcePosition_Type;
TF.SourcePosition_Units=IR.SourcePosition_Units;
TF.EmitterPosition=IR.EmitterPosition;
TF.EmitterPosition_Type=IR.EmitterPosition_Type;
TF.EmitterPosition_Units=IR.EmitterPosition_Units;
TF.ReceiverPosition=IR.ReceiverPosition;
TF.ReceiverPosition_Type=IR.ReceiverPosition_Type;
TF.ReceiverPosition_Units=IR.ReceiverPosition_Units;

TF.Data.Real=zeros(IR.API.M,IR.API.R,IR.API.N+1);
TF.Data.Imag=zeros(IR.API.M,IR.API.R,IR.API.N+1);
for ii=1:IR.API.M
  for jj=1:IR.API.R
   sp=fft(squeeze(IR.Data.IR(ii,jj,:)),2*IR.API.N); % Delay not considered!
   TF.Data.Real(ii,jj,:)=real(sp(1:IR.API.N+1,:));
   TF.Data.Imag(ii,jj,:)=imag(sp(1:IR.API.N+1,:));
  end
end
TF.N=(0:fs/2/IR.API.N:fs/2)';

TF=SOFAupdateDimensions(TF);

%% Convert to an emitter-based representation, TFE
TFE=TF; 
TFE.GLOBAL_SOFAConventions = 'GeneralTF-E';
TFE.GLOBAL_DataType = 'TF-E';
TFE.API.E=TF.API.M;
TFE.API.M=1;
TFE.Data=rmfield(TFE.Data,{'Real','Imag'});
TFE.Data.Real(1,:,:,:)=shiftdim(TF.Data.Real,1); % MRN --> 1RNM --> MRNE with M=1
TFE.API.Dimensions.Data.Real='MRNE';
TFE.Data.Imag(1,:,:,:)=shiftdim(TF.Data.Imag,1);
TFE.API.Dimensions.Data.Imag='MRNE';
TFE.EmitterPosition=TF.SourcePosition;
TFE.EmitterPosition_Type=TF.SourcePosition_Type;
TFE.EmitterPosition_Units=TF.SourcePosition_Units;
TFE.API.Dimensions.EmitterPosition='ECI';
TFE.SourcePosition=[0 0 0];
TFE.API.Dimensions.SourcePosition='IC';

TFE=SOFAupdateDimensions(TFE);


%% Convert to SH
SH=TFE;
SH.GLOBAL_SOFAConventions = 'FreeFieldHRTF';

Lmax=floor(sqrt(size(SH.EmitterPosition,1))-1); % Max SH order
L=2; % actual SH order
[S, SH.API.E]=sph2SH(SH.EmitterPosition(:,1:2), L);

Sinv=pinv(S);
SH.Data.Real=zeros(1, SH.API.R, SH.API.N, SH.API.E);
SH.Data.Imag=zeros(1, SH.API.R, SH.API.N, SH.API.E);
for ii=1:TFE.API.R
  for jj=1:TFE.API.N
   SH.Data.Real(1,ii,jj,:)=Sinv*squeeze(TFE.Data.Real(1,ii,jj,:));
   SH.Data.Imag(1,ii,jj,:)=Sinv*squeeze(TFE.Data.Imag(1,ii,jj,:));
  end
end

SH.EmitterPosition=mean(SH.EmitterPosition);
SH.EmitterPosition_Type='Spherical Harmonics';

SH = SOFAupdateDimensions(SH);


%% interpolate for the horizontal and median planes to SimpleFreeFieldHRTF (TF)
TFint=TF;
ele=[-90:0.5:90 89:-.5:-90 zeros(1,length(1:0.5:355))]';
azi=[zeros(length(-90:.5:90),1); 180*ones(length(89:-.5:-90),1); (1:0.5:355)'];
radius=1.2*ones(size(ele));
TFint.SourcePosition=[azi ele radius];
Sint = sph2SH(TFint.SourcePosition(:,1:2), sqrt(SH.API.E)-1);
TFint.API.M=size(Sint,1);
TFint.Data.Real=zeros(TFint.API.M,2,TFint.API.N);
TFint.Data.Imag=zeros(TFint.API.M,2,TFint.API.N);
for ii=1:TFint.API.R
  for jj=1:TFint.API.N
    TFint.Data.Real(:,ii,jj)=Sint*squeeze(SH.Data.Real(1,ii,jj,:));
    TFint.Data.Imag(:,ii,jj)=Sint*squeeze(SH.Data.Imag(1,ii,jj,:));
  end
end

TFint=SOFAupdateDimensions(TFint);

%% compare
figure;
SOFAplotHRTF(TFint,'magmedian'); title('SimpleFreeFieldHRTF (TF): Interpolated');
figure;
SOFAplotHRTF(IR,'magmedian'); title('SimpleFreeFieldHRIR (FIR) for reference');

% SOFAplotHRTF(TFint,'etchorizontal'); title('SimpleFreeFieldHRTF (TF): Interpolated');

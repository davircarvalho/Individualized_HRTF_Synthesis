function Obj = sofaIRfilter(Obj, fmin, fmax)
fs = Obj.Data.SamplingRate;
IR = shiftdim(Obj.Data.IR,2);

%% Método 1
% % d  = fdesign.bandpass('N,F3dB1,F3dB2', 6, fmin, fmax,fs);  %%% Especifica filtro
% % Hd = design(d, 'butter');   %%% Especifica filtro
% % IR_filtered = filter(Hd, IR, 1);  %% Sinal filtrado

hpFilt  = fdesign.highpass('N,F3dB',4, fmin, fs);
lpFilt  = fdesign.lowpass	('N,F3dB',6, fmax, fs);
Hdhigh  = design(hpFilt, 'butter');   %%% Especifica filtro
Hdlow   = design(lpFilt, 'butter');   %%% Especifica filtro

IR_hp    = filter(Hdhigh, IR, 1);  %% Sinal filtrado
IR_hp_lp = filter(Hdlow, IR_hp, 1);  %% Sinal filtrado


%% Método 2 
% hpFilt = designfilt('highpassiir','FilterOrder',4, ...
%          'PassbandFrequency',fmin,'PassbandRipple',0.2, ...
%          'SampleRate',fs);
% 
% lpFilt = designfilt('lowpassiir','FilterOrder',6, ...
%          'PassbandFrequency',fmax,'PassbandRipple',0.2, ...
%          'SampleRate',fs);
% 
% % fvtool(lpFilt);
% 
% IR_hp    = filtfilt(hpFilt, IR);  %% high pass
% IR_hp_lp = filtfilt(lpFilt, IR_hp);  %% low pass

%% Output
Obj.Data.IR = shiftdim(IR_hp_lp, 1); %% Devolve ao objeto 
end

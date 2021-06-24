function Obj = sofaIRfilter(Obj, fmin, fmax)
fs = Obj.Data.SamplingRate;
IR = shiftdim(Obj.Data.IR,2);

%% Método 1
hpFilt  = fdesign.highpass('N,F3dB',4, fmin, fs);
lpFilt  = fdesign.lowpass('N,F3dB',6, fmax, fs);
Hdhigh  = design(hpFilt, 'butter');   %%% Especifica filtro
Hdlow   = design(lpFilt, 'butter');   %%% Especifica filtro

% IR_hp    = filter(Hdhigh, IR, 1);  %% Sinal filtrado
% IR_hp_lp = filter(Hdlow, IR_hp, 1);  %% Sinal filtrado

% Zero-phase
IR_hp    = filtfilt(Hdhigh.sosMatrix, Hdhigh.ScaleValues,  IR);  %% high pass
IR_hp_lp = filtfilt(Hdlow.sosMatrix, Hdlow.ScaleValues, IR_hp);  %% low pass

%% Output
Obj.Data.IR = shiftdim(IR_hp_lp, 1); %% Devolve ao objeto 
end


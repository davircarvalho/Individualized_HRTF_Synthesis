function Obj = sofaAddLowFrequency(Obj)
fs = Obj.Data.SamplingRate; %sampling frequency
N = size(Obj.Data.IR, 3); %original length of HRTFs
f = (0:N-1)/N*fs; %frequency vector
indf = find(f>=100,1,'first');
IR = shiftdim(Obj.Data.IR,2);
% sweep to expand lower frequencies
flattemp = ita_generate_sweep('samplingRate', fs,...
                              'freqRange', ([5 f(indf)]),...
                              'mode', 'exp');
flat = itaAudio(ifft(flattemp.freqData, N), fs, 'time');
for k = 1:size(IR, 2)
    for l = 1:size(IR, 3)
        IRita = itaAudio(IR(:,k,l),fs,'time');
        IR_outITA = ita_merge(flat,IRita);
        IR_outITA = ita_filter_peak(IR_outITA, 'Q', 0.5, 'fc', 40, 'gain', 12);
        Obj.Data.IR(k,l,:) = IR_outITA.time(:,2);
    end
end
end
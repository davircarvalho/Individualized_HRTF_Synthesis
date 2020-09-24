function itd = sofaGetITD(Obj, varargin)
% Determinação do ITD para HRTFs SOFA
% Davi R. Carvalho @UFSM - Engenharia Acustica - Setembro/2020


%% Parse inputs
% Opção de saída em samples ou tempo
defaultMode = 'samples';
validOutputs = {'samples','time'};
checkOutMode = @(x) any(validatestring(x, validOutputs));

p = inputParser;
addRequired(p,'Obj',@isstruct);
addOptional(p,'outputMode', defaultMode,checkOutMode)

parse(p, Obj, varargin{:})


    %% Comum 
    itd = zeros(length(Obj.SourcePosition), 1);
    IR = shiftdim(Obj.Data.IR, 2);
    fs = Obj.Data.SamplingRate;
          
    %% Metodo 1 : Threshold (usando ita)
    % a partir da diferença entre picos de RI (ISO 3382 A.3.4.).
    A = itaAudio; A.samplingRate = fs;
    B=A;
    for k = 1:size(IR, 2)
        A.time = IR(:,k,1); % L
        B.time = IR(:,k,2); % R    
%         ponto de chegada 
        OnSetL = ita_start_IR(A, 'correlation', false, 'threshold', 10);
        OnSetR = ita_start_IR(B, 'correlation', false, 'threshold', 10);
%         diferenca interaural
        itd(k) = abs(OnSetL - OnSetR);
    end


    %% Metodo 2: phase delay
%     hInput = fft(IR);
%     N = size(hInput, 1);
%     wfreq = 2*pi*(linspace(0, (fs-fs/N), N)).';
%     wfreq(1) = [];  % remove DC
% 
%     for k = 1:size(IR, 2)
%         gama = unwrap(angle(squeeze(hInput(:,k,:))));
%               
%         gama(1,:) = []; % remove DC
%         
%         phasediff = gama(:,1) - gama(:,2); 
% 
%         TOA = mean(phasediff./w); % tempo de chegada (cada canal)
%         
%         itd(k) =  abs(TOA*fs);
%     end

% % % %%% experimental interp to smooth (cuidado)
% % % figure()
% % % plot(itd); hold on 
% % % [peaksu, locsu] = findpeaks(itd, 'MinPeakDistance', 20);
% % % [peaksd, locsd] = findpeaks(-itd, 'MinPeakDistance', 30);
% % % peaks = [itd(1); peaksu; peaksd; itd(end)];
% % % loc = [1; locsu; locsd; length(itd)];
% % % 
% % % plot(loc, peaks, 'x');
% % % 
% % % 
% % % Xq = 1:length(itd);
% % % Vq = abs(interp1(loc,peaks,Xq,'spline'));
% % % plot(Vq, 'linewidth', 1.3)
% % % 
% % % itd = Vq;



%% Método 3  (correlação cruzada com fase mínima)
%     % To find the minimum phase HRIR using Hilbert transform
%     % The minimum phase response and magnitude spectrum are connected
%     % via the Hilbert transform
%     % phaseRespMin = imaginary(HilbertTransform(-ln(magnitude of hSignal)))
%     hSignal = fft(IR);
% 
%     hMag = abs(hSignal); % magnitude        
%     hLnMag = log(hMag);            
%     hHilbert = hilbert(-hLnMag);            
%     hMinPhase = (imag(hHilbert)); % phase
% 
%     hMin = ifft((hMag.*exp(1i*hMinPhase)));
%     hMin = real(hMin);
% 
% %     Obj.Data.IR = shiftdim(hMin, 1);
%     % To find the delay tau
%     % For the estimation of the absolute delay tau , using the spectra
%     % of an impulse response (H[]) and its minimum phase version Hmin[] for
%     % calculating the cross-correlation
%     hMinFFT = fft(hMin);
% 
%     [~, toa] = max(ifft(hSignal.*conj(hMinFFT)));      
%     toa = squeeze(toa);
%     itd = abs(toa(:,1) - toa(:,2));



%% Output Units
switch p.Results.outputMode
    case 'time'
        itd = (itd./Obj.Data.SamplingRate);
end

end

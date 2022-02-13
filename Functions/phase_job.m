function [IR_L, IR_R] = phase_job(hL, hR, itd, pos, fs, offset)
%%% input %%%%
% hL, hR: magnitudes para orelha esquerda e direita (log)
% itd: diferença iteraural de tempo
% azi: azimute referente ao determinado itd 
% pos: azimute (0° -> 360° [SOFA1.0])
%offset: atraso para ambas as RIs
if nargin < 6
    offset = 20.4;
end

%% Log 2 Linear
hL = 10.^(hL./20);
hR = 10.^(hR./20);

hL = [hL; flipud(hL)]; % espelha magnitude
hR = [hR; flipud(hR)];
N = length(hL);

                %%% PHASE RECONSTRUCTION %%%             
%% Minimum Phase 
% Calculo da fase
phi_minL = imag(hilbert(-(log(abs(hL)+eps))));
phi_minR = imag(hilbert(-(log(abs(hR)+eps)))); 

% HRTF complexa
HminL = abs(hL).*exp(1i*phi_minL);
HminR = abs(hR).*exp(1i*phi_minR);

%% Excess phase (ITD) 
if pos(1) >= 180  || (pos(2) < 0 && pos(1) >= 180)
    leftDelay  = (itd+offset);
    rightDelay = (offset); 
else
    leftDelay  = (offset);
    rightDelay = (itd+offset);
end

%% Excess phase (ITD)
HRTF_L = (HminL.'); 
HRTF_R = (HminR.'); 

% Aplicar ITD ainda no dominio da frequencia
% k = 0:N-1;
% LdelayConstant = exp(-1i*2*pi*k*leftDelay / N);
% RdelayConstant = exp(-1i*2*pi*k*rightDelay / N);
% HRTF_L = LdelayConstant .* HRTF_L;
% HRTF_R = RdelayConstant .* HRTF_R;


%%% Back to time domain
IR_L = real(ifft(HRTF_L, N))';
IR_R = real(ifft(HRTF_R, N))';

% Aplicar ITDs inteiros (deslocamento por sample)
% IR_L = circshift(IR_L,leftDelay);
% IR_R = circshift(IR_R,rightDelay);

% Aplicar ITD fracionarios ("interpolar samples intermediarios")
vfd = dsp.VariableFractionalDelay();
IR_L = vfd(IR_L, leftDelay);
IR_R = vfd(IR_R, rightDelay);


end
function [IR_L, IR_R] = phase_job(hL, hR, itd, pos, offset)
%%% input %%%%
% hL, hR: magnitudes para orelha esquerda e direita (log)
% itd: diferença iteraural de tempo
% azi: azimute referente ao determinado itd 
% pos: azimute (0° -> 360° [SOFA1.0])
%offset: atraso para ambas as RIs
if nargin < 5
    offset = 0;
end
offset = offset+itd;

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
if pos(1) >= 180 || (pos(2) < 0 && pos(1) >= 180)
    leftDelay  = round(itd+offset);
    rightDelay = round(offset); 
else
    leftDelay  = round(offset);
    rightDelay = round(itd+offset);
end

%% Excess phase (ITD)
HRTF_L = (HminL.'); 
HRTF_R = (HminR.'); 

%%% Back to time domain
IR_L = real(ifft(HRTF_L, N));
IR_R = real(ifft(HRTF_R, N));

IR_L = circshift(IR_L,leftDelay);
IR_R = circshift(IR_R,rightDelay);

end
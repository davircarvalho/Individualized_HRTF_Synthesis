function [IR_L, IR_R] = phase_job(hL, hR, itd, pos, offset)
%%% input %%%%
% hL, hR: magnitudes para orelha esquerda e direita (log)
% itd: diferença iteraural de tempo
% azi: azimute referente ao determinado itd 
% pos: azimute (0° -> 360° [SOFA1.0])
%offset: atraso para ambas as RIs

%% Log 2 Linear
hL = 10.^(hL./20);
hR = 10.^(hR./20);
N = length(hL);
                %%% PHASE RECONSTRUCTION %%%             
%% Minimum Phase 
% Calculo da fase
phi_minL = imag(hilbert(-log(abs(hL))));
phi_minR = imag(hilbert(-log(abs(hR)))); 

% HRTF complexa
HminL = hL.*exp(1i*phi_minL);
HminR = hR.*exp(1i*phi_minR);

%% Excess phase (ITD) 
if pos(1) >= 180 || (pos(2) < 0 && pos(1) >= 180)
    leftDelay  = itd+offset;
    rightDelay = offset; 
else
    leftDelay  = offset;
    rightDelay = itd+offset;
end

%% Excess phase (ITD)
k = 0:N-1;
LdelayConstant = exp(-1i*pi*k*leftDelay / N);
RdelayConstant = exp(-1i*pi*k*rightDelay / N);

HRTF_L = (HminL.') .* LdelayConstant ;
HRTF_R = (HminR.') .* RdelayConstant ;

%%% Back to time domain
IR_L = real(ifft(HRTF_L, N*2, 'nonsymmetric'));
IR_R = real(ifft(HRTF_R, N*2, 'nonsymmetric'));
end
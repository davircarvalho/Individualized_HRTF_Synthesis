function Obj = hrtf2DTF(Obj, no_samples)
% Transformar HRTF para DTFs, método Middlebrooks(1999)
% Davi R. Carvalho - Janeiro/2020

%% HRTF EM LOG SCALE [HRTF_log]
% Aplica fft e log scale
IR = shiftdim(Obj.Data.IR,2);
[~, no_directions, no_channels] = size(IR);
for k = 1: no_channels
    for l = 1:no_directions
        hrtf_log(:,l,k) = log(abs(fft(IR(:,l,k), no_samples)./no_samples)); 
    end
end
%% CTF
CTF = mean(hrtf_log,2);  


%% DTF
DTF = ((hrtf_log) + (CTF));  

% figure()
% semilogx((CTF(:,1))); hold on
% semilogx((hrtf_log(:,609,1)));
% semilogx((DTF(:,609,1)));
% legend('ctf', 'hrtf', 'dtf')

IR_min = zeros(size(hrtf_log));
for l = 1:no_directions
    %%% Fase mínima %%%
    N = 1:no_samples/2;
    hL = exp(DTF(N, l, 1));
    hR = exp(DTF(N, l, 2));
    phi_minL = imag(hilbert(-log(abs(hL))));
    phi_minR = imag(hilbert(-log(abs(hR)))); 
    % HRTF complexa (mag+phase)
    HminL = hL.*exp(1j*deg2rad(phi_minL));
    HminR = hR.*exp(1j*deg2rad(phi_minR));
    % Back to time domain
    IR_min(:,l,1) = real(ifft(HminL, no_samples, 'symmetric'));
    IR_min(:,l,2) = real(ifft(HminR, no_samples, 'symmetric'));
end
%% OUTPUT
Obj.Data.IR = shiftdim(IR_min,1);
end
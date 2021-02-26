clc; clear all;
% test distnce variation influencce in air absorption coef
N = 200;
fs = 44100;
f = linspace(0, fs-fs/N,N)';
for k = 1:N
    [~, alpha_iso(k,1), ~, ~]=air_absorption(f(k));
end
dist= 5;
alpha = (dist*alpha_iso);
alpha(:,2) = alpha(:,1);
% figure()
% plot(f(1:N/2), alpha, 'r')
% xlim([0 f(N/2)])

%% Aplicar em HRTF 
Obj = SOFAload('..\individuo_141.sofa');
[itd, Obj] = SOFAgetITD(Obj);
ir = shiftdim(Obj.Data.IR, 2);
IR = [ir(:,1,1), ir(:,1,2)];

hrtf = 20*log10(abs(fft(IR)));
y = hrtf - alpha;

% figure()
% plot(f, hrtf); hold on
% plot(f, y)
% legend('Original', 'Corrigido', 'location', 'best')
% title(['Correção da absorção do ar para fonte a ' num2str(dist) 'm'])
    
% plot(y(1:N/2,:))
y_min = (get_min_phase(y(1:N/2,:), 'log', 'nonsymmetric'));
fymin = db(abs((y_min)));
plot( real(ifft(y_min)));

% xlim([0, f(N/2)])


%%
function y = fast_conv(x1, x2)
            nfft = length(x1) + length(x2) - 1;
            y = real(ifft(fft(x1, nfft).*fft(x2, nfft)));
end

% LSD - log spectral distance
%%% REFERECIA: Deep Neural Network Based HRTF Personalization Using
%Anthropometric Measurements
function lsd = LSD(hrir_msrd, hrir_sim, azi)
no_samples = length(hrir_msrd);
k1 = 1;
k2 = no_samples/2;

hrtf_msrd = abs(fft(hrir_msrd, no_samples));
hrtf_sim  = abs(fft(hrir_sim, no_samples));
hrtf_msrd = (hrtf_msrd(1:k2));
hrtf_sim  = (hrtf_sim(1:k2));

%%%% LSD %%%%
soma = sum((20*log10(hrtf_msrd./hrtf_sim)).^2);
lsd = abs(sqrt(soma/(k2 - k1 + 1)));

% %PLOT
% figure()
% semilogx( (abs(hrtf_msrd))); hold on 
% semilogx( (abs(hrtf_sim))); hold off
% legend('msrd', 'sim', 'location', 'best')
end
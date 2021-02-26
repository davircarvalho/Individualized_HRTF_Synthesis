clear all; clc; close all
% Testando se aumentar profuncidade dos vales na simulaçao melhoram
% distorção espectral ou nao
%%
Obj = SOFAload('individuo_141.sofa');
IR = shiftdim(Obj.Data.IR,2);
mag = (abs(fft(IR)));
N = Obj.API.N;
fs = Obj.Data.SamplingRate;
freq = linspace(0,fs-fs/N, N);

% ver a magnitude
% figure()
% plot(freq(1:N/2), mag(1:N/2,1,1)); hold on

% pegar os vales
% indices a partir de 7kHZ
idx_fmin = dsearchn(freq', 7000);
idx_fmax = dsearchn(freq', 17000);

for k = 1:size(mag, 2)
    for l = 1:size(mag, 3)
        mg = mag(idx_fmin:idx_fmax, k, l);   
               
        [peak, loc] = findpeaks(-mg, 'MinPeakDistance', 5);
        loc_norm = idx_fmin+loc.';
        % aumentar amplitude dos vales
        mag(loc_norm, k, l) = abs(mag(loc_norm,k,l)) * (1/2);
        % interporlar valores ao redor
        X = 1:size(mag(1:N/2,k,l));
        Xq = X;
        V = mag(1:N/2,k,l);
        rm = [];
        for m = 1:length(loc_norm)
            temp = [loc_norm(m)-3: loc_norm(m)-1,...
                 loc_norm(m)+1: loc_norm(m)+3];
            rm = [rm temp];
        end
        X(rm) = [];
        V(rm) = [];     
        Vq = interp1(X,V,Xq, 'spline').';
        
        % pegar fase
        IR_fix(:, k, l) = real(ifft(get_min_phase(Vq,'linear','nonsymmetric')));
%         plot(freq(idx_f+loc), -peak, 'x');
%         plot(freq(1:N/2), mag(1:N/2,1,1), 'linewidth', 1.3)
%         pause(0)
    end
end

Obj2 = Obj;
Obj2.Data.IR = shiftdim(IR_fix, 1);

% PLOT plano horizontal
figure()
type = 'MagHorizontal';
SOFAplotHRTF(Obj, type);
title('DTFs medidas')
xlabel('Frequência [Hz]')
ylabel('Azimute')
axis([100 19000, -175 180])
yticks([-90, 0, 90])
yticklabels({'270°', '0°', '90°'})
set(gca,'FontSize',11)
% export_fig([pwd, '\Images\hrtf_horizontal_med'], '-pdf', '-transparent');

hFigure = figure();
SOFAplotHRTF(Obj2, type);
title('DTFs simuladas')
xlabel('Frequência [Hz]')
ylabel('Azimute')
axis([100 19000, -175 180])
yticks([-90, 0, 90])
yticklabels({'270°', '0°', '90°'})
set(gca,'FontSize',11)

filename = [pwd, '\Images\hrtf_horizontal_sim.pdf'];
% exportgraphics(hFigure,filename,'BackgroundColor','none','ContentType','vector')

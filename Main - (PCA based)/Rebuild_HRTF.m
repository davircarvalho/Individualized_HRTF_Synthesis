clear all; close all; clc
% DAVI ROCHA CARVALHO; ENG. ACUSTICA - UFSM; JUNHO/2019
% Simulação de HRTF individualizada, a partir da RNA treinada e uso de dados
% antropométricos
addpath(genpath([pwd, '\..\Functions'])); 
addpath(genpath([pwd, '\..\EAC-Toolbox\exportfig']));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Options 
% Defina quais datasets usar: {'cipic', 'ari', 'ita', '3d3a', 'riec'}
Datasets = {'cipic', 'ari', 'ita', '3d3a'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load networks
net_file = 'DADOS_TREINAMENTO/net_treinada';
target_file = 'DADOS_TREINAMENTO/target_pca';
if any(strcmp('cipic', Datasets)) 
    net_file = append(net_file, '_CIPIC');
    target_file = append(target_file, '_CIPIC');
end
if any(strcmp({'ari'}, Datasets)) 
    net_file = append(net_file, '_ARI');
    target_file = append(target_file, '_ARI');
end
if any(strcmp({'ita'}, Datasets)) 
    net_file = append(net_file, '_ITA');
    target_file = append(target_file, '_ITA');
end
if any(strcmp({'3d3a'}, Datasets)) 
    net_file = append(net_file, '_3D3A');
    target_file = append(target_file, '_3D3A');
end
if any(strcmp({'riec'}, Datasets)) 
    net_file = append(net_file, '_RIEC');
    target_file = append(target_file, '_RIEC');
end
if any(strcmp({'tub_meas'}, Datasets)) 
    net_file = append(net_file, '_TUBMEAS');
    target_file = append(target_file, '_TUBMEAS');
end
if any(strcmp({'tub_sim'}, Datasets)) 
    net_file = append(net_file, '_TUBSIM');
    target_file = append(target_file, '_TUBSIM');
end

load(net_file)
load(target_file)
fs = 44100;

disp('Network carregada!')

%% Definir INPUTs PRO TESTE
% indivíduos para VALIDAÇÃO: 1('003') .
subj = 1;
% nome do individuo dentro do banco de dados
subject = '003';

% separando os 8 parametros necessários para a ann
anthro = load('anthro_CIPIC.mat');
% d1, d3 : d6
d1 = anthro.D(subj, [1,2,3,5,7,8]);      % L
% d1 = anthro.D(subj, [1,3,5,6,7,8]);      % L
d2 = anthro.D(subj, [9,10,11,13,15,16]);    % R
% d2 = anthro.D(subj, [9,11,13,14,15,16]);    % R

% x1, x3, x12
x = anthro.X(subj, [1,3]);

% unir as duas matrizes
InptMtx(:,:,1) = [x, d1]';  % L 
InptMtx(:,:,2) = [x, d2]';  % R

[no_samples, no_PC, no_directions, no_channels] = size(PCWs);

%% Preprocessamento input
% O desvio e a média utilizados são oriundos do dataset de treinamento
% [xinpt, ~, ~] = nn_preprocess(InptMtx, 'apply', 'sig', sig1, 'mu', mu1);
xinpt = InptMtx;
%% Simulação do modelo
disp('...')
disp('Simulação iniciada')
% computar valores de saida para dada entrada na ann
tic
wait = waitbar(0,'Processando Novas Componentes Principais');
cont = 0;
for n = 1:no_channels
    for i = 1:no_directions
        result = net_pca{i, n}(xinpt(:,:, n));
        DTF_sim(:, i, n) = (PCWs(:, :, i, n) * result) + med_vec2(:, i, n); 
        
        cont = cont+1;
        waitbar((cont)/(no_directions*no_channels), wait);        
    end
end
toc
close(wait)
disp('Simulação concluida!')
disp('...')


%% RECONSTRUÇÃO DE FASE: (fase mínima + ITD)
% CÁlCULO DO ITD
new_itd = itd_synthesis(x(1), x(2), out_pos, fs, 'adapt');

% Aplicando fase minima e de excesso 
for k = 1:no_directions
    itd = new_itd(k);
    offset = 20;  
    [IR_minL, IR_minR] = phase_job(DTF_sim(:, k, 1), DTF_sim(:, k, 2), ...
                                   itd, out_pos(k,:), offset); 
    %save no formato CIPIC
    hrir_final(:, k, 1) = IR_minL;
    hrir_final(:, k, 2) = IR_minR;
end


%% Export para SOFA
% SIMULATED DATA
Obj_sim = SOFAgetConventions('SimpleFreeFieldHRIR');
Obj_sim.Data.IR = shiftdim(hrir_final, 1);
Obj_sim.Data.SamplingRate = fs;
Obj_sim.SourcePosition = out_pos;
Obj_sim = SOFAupdateDimensions(Obj_sim);


% plot(new_itd); hold on; plot(sofaGetITD(Obj_sim));
% Preprocessamento do database para comparação
pathcipic = dir('B:\Documentos\#3 - TCC\EAC-TCC-Davi\HRTF-datasets\Banco CIPIC\SOFA\*.sofa');
CIPIC = SOFAload([pathcipic(subj).folder '\' pathcipic(subj).name], 'nochecks');
% HRTF -> DTF
Obj_med = sofaFit2Grid(CIPIC, out_pos);   

fmin = 250; fmax = 18000;
Obj_med = sofaIRfilter(Obj_med, fmin, fmax);
Obj_med = SOFAhrtf2dtf(Obj_med); 


itdori = sofaGetITD(Obj_med);
res = mean(abs(itdori - sofaGetITD(Obj_sim)))

%% % SAVE %%%
file = 'DADOS_TREINAMENTO/rebuilt_data';
if any(strcmp('cipic', Datasets)) 
    file = append(file, '_CIPIC');
end
if any(strcmp({'ari'}, Datasets)) 
    file = append(file, '_ARI');
end
if any(strcmp({'ita'}, Datasets)) 
    file = append(file, '_ITA');
end
if any(strcmp({'3d3a'}, Datasets)) 
    file = append(file, '_3D3A');
end
if any(strcmp({'riec'}, Datasets)) 
    file = append(file, '_RIEC');
end
if any(strcmp({'tub_sim'}, Datasets)) 
    file = append(file, '_TUBSIM');
end
if any(strcmp({'tub_meas'}, Datasets)) 
    file = append(file, '_TUBMEAS');
end
save([file '.mat'], 'Obj_sim', 'Obj_med')
disp('SOFA files saved!')


%%%% save SOFA file %%%%
% file_path = fullfile(['individuo_' subject '.sofa']);
% compression = 0;
% SOFAsave(file_path, Obj_sim, compression);
























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT - HRTF_log para determinada direção º
N = size(Obj_med.Data.IR, 3);  
freq = linspace(0, fs-fs/N, N);
f_1000hz = dsearchn(freq', 500); %posicao da sample em 1000Hz
Obj_med.SourcePosition = round(Obj_med.SourcePosition);

azimute = [0, 30, 90, 120, 220, 240, 270, 330];
for i = 1:length(azimute) % plot loop
%%%% DEFINA A DIREÇÃO DA HRTF %%%
azim = azimute(i); 
elev = 0;

% Get index of measurements with the same directions
pos=find(Obj_med.SourcePosition(:,1)==azim & Obj_med.SourcePosition(:,2)==elev);


%%% Medido ----------------------------------------------------------------
% L
nfft = N;
gL = squeeze(Obj_med.Data.IR(pos, 1, :));
g_L = fft(gL, nfft); 
g_L=g_L./g_L(f_1000hz); %normalizar em 1000Hz
g_L = 20*log10(abs(g_L));

% R 
gR = squeeze(Obj_med.Data.IR(pos, 2, :));
g_R = fft(gR, nfft); 
g_R = g_R./g_R(f_1000hz); %normalizar em 1000Hz
g_R = 20*log10(abs(g_R));

%%% Simulado --------------------------------------------------------------
% L
hrtf_simL = squeeze(Obj_sim.Data.IR(pos, 1, :));
h_simL = fft(hrtf_simL, nfft);
h_simL = h_simL./h_simL(f_1000hz); %normalizar em 1000Hz
h_simL = 20*log10(abs(h_simL));
% R 
hrtf_simR = squeeze(Obj_sim.Data.IR(pos, 2, :));
h_simR = fft(hrtf_simR, nfft); 
h_simR = h_simR./h_simR(f_1000hz); %normalizar em 1000Hz
h_simR = 20*log10(abs(h_simR));



%%%% PLOT %%%%
figure()
semilogx(freq(1:N/2), g_L(1:N/2), 'lineWidth', 1.5,  'Color', 'blue'); hold on
semilogx(freq(1:N/2), g_R(1:N/2), 'lineWidth', 1.5, 'Color', 'red'); 


semilogx(freq(1:N/2), h_simL(1:N/2), '-.r', 'lineWidth', 1.5, 'Color', 'blue'); 
semilogx(freq(1:N/2), h_simR(1:N/2), '-.r', 'lineWidth', 1.5, 'Color', 'red'); hold off

% PT_BR
% legend('Original Esquerda', 'Original Direita', ... 
%        'Proposta Esquerda', 'Proposta Direita', ...
%                'Location', 'best', 'lineWidth', 1.5)
%     
% xlabel('Frequência [Hz]');
% ylabel('Amplitude [dB]');
% EN_US
legend('Original left', 'Original right', ... 
       'Predicted left', 'Predicted right', ...
               'Location', 'best', 'lineWidth', 1.5)
    
xlabel('Frequency [Hz]');
ylabel('Amplitude [dB]');


azim = num2str(out_pos(pos,1));
elee = num2str(out_pos(pos,2));
% title(['DTF Individuo ', subject,', azimute ' azim 'º, elevação ' elee 'º.' ]);
title(['DTF Subject ', subject,', azimuth ' azim 'º, elevation ' elee 'º.' ]);


grid on 
xlim([250 17000])
set(gca, 'YLimSpec', 'Tight');
% axis tight
set(gca,'FontSize',13)


% export_fig([pwd, '\Images\English\hrtf' num2str(azimute(i))], '-pdf', '-transparent');
end
%% PLOT - HRIR [Simulada vs Medida]
figure()

% DEFINA A DIREÇÃO DA RI
azim = 270; 
elev = 0;

% Get index of measurements with the same directions
pos=find(round(Obj_med.SourcePosition(:,1))==azim & round(Obj_med.SourcePosition(:,2))==elev);


%%% ESQUERDA --------------------------------------------------------------
subplot(2,1,1)
% Medida
plot(real(squeeze(Obj_med.Data.IR(pos, 1, :))), 'linewidth', 1.5);hold on
% Proposta
plot(real(squeeze(Obj_sim.Data.IR(pos, 1, :))), '-.r', 'linewidth', 1.3); hold off

xlabel('Amostras')
ylabel('Esquerda')
azim = num2str(out_pos(pos,1));
elee = num2str(out_pos(pos,2));
title(['HRIR Invididuo ', num2str(subject), ', azimute ', azim, ...
                                   '°, elevação ', elee, '°']);
legend('Original', 'Simulada')
set(gca,'FontSize',12)
axis tight

%%% DIREITA ---------------------------------------------------------------
subplot(2,1,2)
plot( real(squeeze(Obj_med.Data.IR(pos, 2, :))), 'linewidth', 1.5); hold on
plot( real(squeeze(Obj_sim.Data.IR(pos, 2, :))) , '-.r','linewidth', 1.3); hold off

xlabel('Amostras')
ylabel('Direita')
legend('Original', 'Simulada')

% ylim([-1 1])
set(gca,'FontSize',12)
axis tight




%% PLOT plano horizontal
figure()
type = 'MagHorizontal';
SOFAplotHRTF(Obj_med, type);
title('DTFs medidas')
xlabel('Frequência [Hz]')
ylabel('Azimute')
axis([100 19000, -175 180])
yticks([-90, 0, 90])
yticklabels({'270°', '0°', '90°'})
set(gca,'FontSize',11)
% export_fig([pwd, '\Images\hrtf_horizontal_med'], '-pdf', '-transparent');

hFigure = figure();
SOFAplotHRTF(Obj_sim, type);
title('DTFs simuladas')
xlabel('Frequência [Hz]')
ylabel('Azimute')
axis([100 19000, -175 180])
yticks([-90, 0, 90])
yticklabels({'270°', '0°', '90°'})
set(gca,'FontSize',11)

filename = [pwd, '\Images\hrtf_horizontal_sim.pdf'];
% exportgraphics(hFigure,filename,'BackgroundColor','none','ContentType','vector')


%% PLOT plano vertical
% figure
% type = 'MagSagittal';
% subplot(211)
% SOFAplotHRTF(Obj_med, type);
% title('DTFs medidas')
% xlabel('Frequência [Hz]')
% ylabel('Azimute')
% % axis([250 19000, -175 180])
% % yticks([-90, 0, 90])
% % yticklabels({'270°', '0°', '90°'})
% set(gca,'FontSize',11)
% 
% subplot(212)
% SOFAplotHRTF(Obj_sim, type);
% title('DTFs simuladas')
% xlabel('Frequência [Hz]')
% ylabel('Azimute')
% % axis([250 19000, -175 180])
% % yticks([-90, 0, 90])
% % yticklabels({'270°', '0°', '90°'})
% set(gca,'FontSize',11)
% export_fig([pwd, '\Images\hrtf_vertical'], '-pdf', '-transparent');



















%% Low Frequeny Extension [LFE] (TO DO)

% fs = 44100;
% Ni = 4096;  
% fmax = 315;
% freq = linspace(0, fs-fs/Ni, Ni);
% LFE_mag = ones(Ni/2, 1); 
% LFE_time = ifft(LFE_mag);
% % LFE minimum phased 
% phi_min = imag(hilbert(-log(abs(LFE_mag))));
% Hmin = LFE_mag.*exp(1i*phi_min);
% LFE_time = real(ifft(Hmin, Ni, 'nonsymmetric'));
% 
% lfe_time = db(abs(fft(LFE_time)));
% % figure()
% % semilogx(freq(2:Ni/2), lfe_time(2:Ni/2)); hold on
% 
% % FILTRAGEM 
% % d  = fdesign.bandpass('N,F3dB1,F3dB2',2, 5, fmax,fs);  %%% Especifica filtro
% d  = fdesign.lowpass('N,F3dB',2, fmax, fs);
% Hdlow = design(d, 'butter');   %%% Especifica filtro
% LFE = filter(Hdlow, LFE_time, 1);  %% Sinal filtrado
% % lfe = db(abs(fft(LFE)));
% % semilogx(freq(2:Ni/2), lfe(2:Ni/2)); hold off
% % axis([0 1000 -10 -5])
% 
% 
% % ajustar nivel da LFE de acordo com low freq do input
% irf = DTF_sim(:,1,1);
% 
% lfe_freq = abs(fft(LFE));
% lfe_freq(1) = lfe_freq(2);
% 
% f_cross = dsearchn(freq', 200); % cross over freq
% 
% novo_nivel = irf(f_cross+40);
% 
% lfe_freq = (lfe_freq./(max(lfe_freq(:)))).*novo_nivel;
% 
% comb = [lfe_freq(1:f_cross); irf(f_cross:end)];
% 
% 
% % interpolate
% X = [1:f_cross-10, f_cross+20:f_cross+50];
% V = comb(X);
% Xq = (f_cross-10 : f_cross+25);
% Vq = interp1(X,V,Xq, 'spline');
% comb(Xq) = Vq; 
% 
% % MAGNITUDE
% figure()
% semilogx(freq(2:Ni/2), db(irf(2:Ni/2)), 'linewidth', 2); hold on 
% semilogx(freq(2:Ni/2), db(lfe_freq(2:Ni/2)), '--r', 'linewidth', 1); 
% semilogx(freq(2:Ni/2), db(comb(2:Ni/2)), 'linewidth', 2);  
% axis([0 20000 -50 -15 ])
% legend('Original', 'LFE', 'Original + LFE', 'Location', 'Southwest')
% 
% 
% %% PHASE
% figure()
% combt = ifft(comb);
% plot(unwrap(phase(ir)), 'linewidth', 2); hold on 
% plot(unwrap(phase(combt)), 'linewidth', 2);  
% % axis([0 20000 -50 -15 ])
% legend('Original', 'Original + LFE', 'Location', 'Southwest')

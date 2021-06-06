clear all; clc; 
% DAVI ROCHA CARVALHO; ENG. ACUSTICA - UFSM; NOVEMBRO/2020
% Analise de componentes principais 
addpath(genpath([pwd, '\..\Functions'])); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Datasets: {'cipic', 'ari', 'ita', '3d3a', 'riec', 'tub_meas', 'tub_sim'}

Datasets = {'cipic', 'ari', 'ita', '3d3a'};
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load networks
DTF_file = [pwd, '\..\DADOS_TREINAMENTO\DTF'];
if any(strcmp('cipic', Datasets)) 
    DTF_file = append(DTF_file, '_CIPIC');
end
if any(strcmp({'ari'}, Datasets)) 
    DTF_file = append(DTF_file, '_ARI');
end
if any(strcmp({'ita'}, Datasets)) 
    DTF_file = append(DTF_file, '_ITA');
end
if any(strcmp({'3d3a'}, Datasets)) 
    DTF_file = append(DTF_file, '_3D3A');
end
if any(strcmp({'riec'}, Datasets)) 
    DTF_file = append(DTF_file, '_RIEC');
end
if any(strcmp({'tub_meas'}, Datasets)) 
    DTF_file = append(DTF_file, '_TUBMEAS');
end
if any(strcmp({'tub_sim'}, Datasets)) 
    DTF_file = append(DTF_file, '_TUBSIM');
end
load(DTF_file)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pre-processsamento
no_samples = size(DTF, 1);
DTF_ok = DTF((1:no_samples/2),:,:,:); % take the single sided spectra
[no_samples, no_subj, no_directions, no_channels] = size(DTF_ok);


%% Principal Component Analysis (PCA)
clear med_vec2 coeffs PCWs explain

no_PC = 12; %número de principais componentes de interesse
for m = 1:no_channels
    for n = 1:no_directions
%         med_vec2(:,n,m) = mean(DTF_ok(:, :, n, m), 2);
        data_mtx = DTF_ok(:, :, n, m);% - med_vec2(:,n,m);
         
        [coeff, score, ~,~,explained,mu] = pca(data_mtx.','NumComponents', no_PC, ...
                                                        'Centered', true, ...
                                                        'Algorith', 'svd');      
        coeffs(:,:, n, m) = coeff;
        % Scores are the representations of X in the principal component space
        PCWs(:,:, n, m) = score;
        explain(:,:,n,m) = explained;
        med_vec2(:,n,m) = mu;
    end
end
disp('PCA calculada!')


%% PLOT DO NUMERO DE PC NECESSÀRIOS 

hFigure = figure();
y = cumsum(mean(explain(:,:,:,1),3));
h = plot(y, 'LineWidth', 2.0);
axis tight
ylabel('Recuperação de variância total [%]')
xlabel('Número de Componentes Principais [~]')
grid on 
idx_tip = dsearchn(y, 90);
datatip(h, idx_tip, y(idx_tip));
set(gca, 'FontSize', 12)

filename = [pwd, '\Images\no_PC.pdf'];
exportgraphics(hFigure,filename,'BackgroundColor','none','ContentType','vector')


%% SAVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path_save = [pwd, '\..\DADOS_TREINAMENTO\target_pca'];
if any(strcmp('cipic', Datasets)) 
    path_save = append(path_save, '_CIPIC');
end
if any(strcmp({'ari'}, Datasets)) 
    path_save = append(path_save, '_ARI');
end
if any(strcmp({'ita'}, Datasets)) 
    path_save = append(path_save, '_ITA');
end
if any(strcmp({'3d3a'}, Datasets)) 
    path_save = append(path_save, '_3D3A');
end
if any(strcmp({'riec'}, Datasets)) 
    path_save = append(path_save, '_RIEC');
end
if any(strcmp('tub_meas', Datasets))
    path_save = append(path_save, '_TUBMEAS');
end
if any(strcmp('tub_sim', Datasets))
    path_save = append(path_save, '_TUBSIM');
end
save(path_save, 'coeffs', 'PCWs', 'med_vec2', 'out_pos')
disp('Dados salvos!')






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot reconstruction
% N = 2*no_samples;  
% freq = linspace(0, fs-fs/N, N);
% subj = 2;
% ch = 2;
% dir = 150;
% recon = PCWs(:,:,dir,ch)*coeffs(subj,:,dir,ch)'+ med_vec2(:,dir,ch);
% 
% figure()
% plot(freq(1:N/2), DTF_ok(:,subj,dir,ch), 'linewidth', 2.5, 'color', [0 0 0]); hold on 
% plot(freq(1, 1:N/2), recon(:,1), 'r','linewidth', 1.5); 
% 
% ylim([-30 5]); xlim([900 2e4])
% axis tight
% xlabel('Frequencia [Hz]'); ylabel('Amplitude [dB]')
% legend('Original', [num2str(no_PC) ' CPs'], 'Location', 'best')
% set(gca, 'FontSize', 13)



%% Avaliação da reconstrução 
clear recon
fmin = 0;
fmax = 20000;
recon = zeros(size(DTF_ok));

for m = 1:no_channels
    for l = 1:no_directions
        for k = 1:no_subj
            recon(:,k,l,m) = coeffs(:,:,l,m)*PCWs(k,:,l,m)'+ med_vec2(:,l,m);            
            SD(:,k,l,m) = spectral_dist(recon(:,k,l,m),DTF_ok(:,k,l,m),fs,fmin,fmax);
        end
    end
end


%% Mapa projeção
ch = 1;
SD_recon_PCA = squeeze(mean(SD,2));
% color_range = [0, 2];
color_range = [0 0.5];

% Simulada 
hFigure = figure('Renderer', 'painters', 'Position', [10 10 600 460]);
scatter(out_pos(:,1), out_pos(:,2), 30, SD_recon_PCA(:,ch), 'filled')

xlabel('Azimuth (°)')
ylabel('Elevation (°)')
title('Case 2')
axis tight
c = colorbar; caxis(color_range); colormap jet
c.Label.String = 'Spectral distortion (dB)';
set(gca,'FontSize',13)
% set(gca,'Color','k')

filename = [pwd, '\Images\MAP_PCAerror.pdf'];
% exportgraphics(hFigure,filename,'BackgroundColor','none','ContentType','vector')


%% Probabilidade 
% hFigure = figure('Renderer', 'painters', 'Position', [10 10 600 460]);
% histogram(SD_recon_PCA(:,ch), 'Normalization','probability', 'NumBins', 13)
% ytix = get(gca, 'YTick');
% set(gca, 'YTick',ytix, 'YTickLabel',ytix*100);
% 
% xlabel('Distorção espectral [dB]')
% ylabel('Probabilidade [%]')
% xticks([0:0.2:2])
% xlim([0, 2])
% set(gca,'FontSize',12)
% filename = [pwd, '\Images\Prob_PCAerror.pdf'];
% exportgraphics(hFigure,filename,'BackgroundColor','none','ContentType','vector')


%% INTERNAL FUNCTIONS %----------------------------------------------------
function sd = spectral_dist(msrd,ref,fs,fmin,fmax)
%% Equals the size of both vectors
N = 2*length(ref);          
%% Create frequency vector to find value in 1000 Hz
f = linspace(0,fs*(N-1)/N,N);
f_500hz = dsearchn(f',500);
fmin = dsearchn(f',fmin);
fmax = dsearchn(f',fmax);
%% Normalizes vectors in frequency domain
msrd = msrd./msrd(f_500hz);
ref = ref./ref(f_500hz);
%% Calculates the Spectral Distortion
sd = sqrt((1/(fmax-fmin+1))*sum((20*log10(abs(msrd(fmin:fmax)./ref(fmin:fmax)))).^2));
end

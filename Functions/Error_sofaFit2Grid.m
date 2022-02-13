clear all; clc
% Avaliação de performance entre modelos de interpolação de posições de HRTFs
% Davi R. Carvalho - Julho/2021

%% LOAD 
local = [pwd '\..\Datasets\'];    
dataset_name = 'HUTUBS';
path = dir([local dataset_name '\pp*_HRIRs_simulated.sofa']);
% % path = dir([local dataset_name '\*simulated.sofa']);
% path = dir([local dataset_name '\*.sofa']);


%% Config
n_subjects = length(path);
dataset = cell(1,n_subjects);
for k = 1:n_subjects
   dataset{k} = SOFAload([path(k).folder, '\',path(k).name], 'nochecks');
   if strcmp(dataset_name, 'AACHEN')
       dataset{k} = ITA2spheric(dataset{k});     
   end
end
if strcmp(dataset_name, 'AACHEN')
    dataset(14) = [];
end
n_subjects = length(dataset);
% n_subjects = 3;
% no_posi = size(dataset{1}.SourcePosition,1); % numero total de posições
idx_pos = find(dataset{1}.SourcePosition(:,2)> -60 & ...
               dataset{1}.SourcePosition(:,2)<  60); % limit elevation range
posi = dataset{1}.SourcePosition(idx_pos, :);
no_posi = length(posi);


%% Main loop
cont=1;
n_desired_pos = [100]; % número de posições "objetivo" (removidas pro teste)
for n = 1:length(n_desired_pos) 
    rng(0) % reset random generator
    idx_temp = randperm(no_posi, n_desired_pos(n)); % Index de posições a serem removidas
    idx = dsearchn(dataset{1}.SourcePosition, posi(idx_temp, :));
    for ii = 1:n_subjects  % numero de individuos sob analise 
        tic
        % HRIRs
        inpt_IR = dataset{ii}.Data.IR; % all HRIR for this subject
        des_IR  = inpt_IR(idx,:,:); % HRIRs objetivo
        inpt_IR(idx,:,:) = [];      % HRIRs de entrada (known positions)

        % Coordinates  
        inpt_pos = dataset{1}.SourcePosition;
        des_pos  = inpt_pos(idx,:); % Posicoes objetivo
        inpt_pos(idx,:)  = [];      % Posicoes de entrada    

        % SOFA sem posições objetivo (ENTRADA)
        Obj_inpt = generate_sofa(inpt_pos, inpt_IR, dataset{ii}.Data.SamplingRate);

        % SOFA apenas com posições objetivo (OBJETIVO)
        REAL{n} = generate_sofa(des_pos, des_IR, dataset{ii}.Data.SamplingRate);

        % error
        [ADPT{n}, sd(n).adpt(:,ii), ITD_error(n).adpt(:,ii), ILD_error(n).adpt(:,ii)] = ...
                                        EvaluateDistortions(Obj_inpt, REAL{n}, 'adapt');
        [VBAP{n}, sd(n).vbap(:,ii), ITD_error(n).vbap(:,ii), ILD_error(n).vbap(:,ii)] = ...
                                        EvaluateDistortions(Obj_inpt, REAL{n}, 'vbap');
        [BLIN{n}, sd(n).blin(:,ii), ITD_error(n).blin(:,ii), ILD_error(n).blin(:,ii)] = ...
                                        EvaluateDistortions(Obj_inpt, REAL{n}, 'bilinear');
        [HSPH{n}, sd(n).hsph(:,ii), ITD_error(n).hsph(:,ii), ILD_error(n).hsph(:,ii)] = ...
                                        EvaluateDistortions(Obj_inpt, REAL{n}, 'spherical_harmonics');
        
        clc; disp(['Iteração: ' num2str(n) ' de ' num2str(length(n_desired_pos)), ....
                   '    Subject: ' num2str(ii), ' de ' num2str(n_subjects)]); 
 
        ETA = toc*(length(n_desired_pos)*n_subjects-cont)/60;
        disp(['ETA: ' num2str(ETA), ' min' ]) 
        cont = cont+1;
        
        
%         plot_mag_map(REAL{n}, HSPH{n})
    end
    %save para plot
    plt(n).des  = des_pos;
    plt(n).inpt = inpt_pos;
end

%% SAVE THEM ALL
save([pwd '\..\DADOS_TREINAMENTO\workspace_sofaFit2Grid_' dataset_name '.mat'])
% clear all; clc
% dataset_name = 'HUTUBS';
% load([pwd '\..\DADOS_TREINAMENTO\workspace_sofaFit2Grid_' dataset_name '.mat'])




    





















%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simplicação do erro para comparação (ITD, ILD)
for k = 1:n
    media_itd(k).adpt = mean(mean([ITD_error(k).adpt], 2));
    media_itd(k).vbap = mean(mean([ITD_error(k).vbap], 2));
    media_itd(k).blin = mean(nanmean([ITD_error(k).blin], 2));
    media_itd(k).hsph = mean(mean([ITD_error(k).hsph], 2));
    
    std_itd(k).adpt = std(mean([ITD_error(k).adpt], 2));
    std_itd(k).vbap = std(mean([ITD_error(k).vbap], 2));
    std_itd(k).blin = std(nanmean([ITD_error(k).blin], 2));
    std_itd(k).hsph = std(mean([ITD_error(k).hsph], 2));
    
    media_ild(k).adpt = mean(mean([ILD_error(k).adpt], 2));
    media_ild(k).vbap = mean(mean([ILD_error(k).vbap], 2));
    media_ild(k).blin = mean(nanmean([ILD_error(k).blin], 2));
    media_ild(k).hsph = mean(mean([ILD_error(k).hsph], 2));
    
    std_ild(k).adpt = std(mean([ILD_error(k).adpt], 2));
    std_ild(k).vbap = std(mean([ILD_error(k).vbap], 2));
    std_ild(k).blin = std(nanmean([ILD_error(k).blin], 2));
    std_ild(k).hsph = std(mean([ILD_error(k).hsph], 2));
end

% figure()
% plot(1:n, [media_itd.adpt]*1e6, 'linewidth', 1.7);hold on
% plot(1:n, [media_itd.vbap]*1e6, 'linewidth', 1.7);
% plot(1:n, [media_itd.blin]*1e6, 'linewidth', 1.7);
% plot(1:n, [media_itd.hsph]*1e6, 'linewidth', 1.7);
% xticks(1:n);
% title('ITD error')
% xticklabels(n_desired_pos)
% xlabel('Estimated positions')
% ylabel('Erro medio (\mus)')
% legend('Nearest', 'VBAP', 'Bilinear',  'Spherical harmonics','location', 'best')
% 

% figure()
% plot([media_ild.adpt], 'linewidth', 1.7);hold on
% plot([media_ild.vbap], 'linewidth', 1.7);
% plot([media_ild.blin], 'linewidth', 1.7);
% plot([media_ild.hsph], 'linewidth', 1.7);
% xticks(1:n);
% title('ild')
% legend('adapt', 'vbap', 'blin',  'hsph','location', 'best')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simplicação do erro para comparação (LSD)
for ks = 1:n
    % média do erro de todas as posições para cada individuo
    pos_sd(ks).adpt = nanmean(sd(ks).adpt);
    pos_sd(ks).vbap = nanmean(sd(ks).vbap);
    pos_sd(ks).blin = nanmean(sd(ks).blin);
    pos_sd(ks).hsph = nanmean(sd(ks).hsph);

    % média de todas as posições e todos os individuos
    media_sd(ks).adpt = mean(pos_sd(ks).adpt);
    stdev(ks).adpt    = std(pos_sd(ks).adpt); %std entre media de cada indivíduo

    media_sd(ks).vbap = mean(pos_sd(ks).vbap);
    stdev(ks).vbap    = std(pos_sd(ks).vbap);

    media_sd(ks).blin = mean(pos_sd(ks).blin);
    stdev(ks).blin    = std(pos_sd(ks).blin);
    
    media_sd(ks).hsph = mean(pos_sd(ks).hsph);
    stdev(ks).hsph    = std(pos_sd(ks).hsph);
end


%% Media e desvio da SD para todas as posições e todos os individuos 
% padrao entre individuos 
% hFigure = figure();
% set(0,'DefaultLineLineWidth',1.5)
% 
% fig1 = shadedErrorBar(1:n, [media_sd.adpt], [stdev.adpt],{'Color',colors(0),'LineWidth', 1.7},1,0.2,'lin', 0);hold on
% fig2 = shadedErrorBar(1:n, [media_sd.vbap], [stdev.vbap],{'Color',colors(1),'LineWidth', 1.7},1,0.2,'lin', 0);
% fig3 = shadedErrorBar(1:n, [media_sd.blin], [stdev.blin],{'Color',colors(2),'LineWidth', 1.7},1,0.2,'lin', 0);
% fig5 = shadedErrorBar(1:n, [media_sd.hsph], [stdev.hsph],{'Color',colors(6),'LineWidth', 1.7},1,0.2,'lin', 0);hold off
% 
% %metadata
% h = get(gca,'Children');
% legendstr = {'', 'Nearest', '', 'VBAP', '', 'Bilinear', '', 'Spherical harmonics'};
% legend(h([1 3 5 7 ]), legendstr{[ 8 6 4 2]}, 'Location', 'Best')
% xticks(1:n)
% xticklabels(n_desired_pos)
% title('')
% xlabel('Estimated positions')
% ylabel('Spectral distortion (dB)')
% set(gca,'FontSize',12)
% ylim([-1 10])
% filename = [pwd, '\Images\ShadedError_sofaFit2Grid.pdf'];
% % exportgraphics(hFigure,filename,'BackgroundColor','none','ContentType','vector')


%% Mapa de ERRO por posição 
lim_colorbar = [1 8.5];
n_idx = 1; % 
plt_pos = plt(n_idx).des;
plt_inp = plt(n_idx).inpt;

% ADPT
hFigure = figure();
scatter(plt_pos(:,1), plt_pos(:,2), 35, nanmean(sd(n_idx).adpt,2), 'filled'); 
title('Nearest')
xlabel('Azimuth (°)')
ylabel('Elevation (°)')
axis tight
c = colorbar; caxis(lim_colorbar); colormap jet
c.Label.String = 'Spectral distortion (dB)';
set(gca,'FontSize',12)
% set(gca,'Color','k');
filename = [pwd, '\Images\MAP_ADAPT_sofaFit2Grid.pdf'];
% exportgraphics(hFigure,filename,'BackgroundColor','none','ContentType','vector')
 
%close  VBAP
hFigure = figure();
scatter(plt_pos(:,1), plt_pos(:,2), 35, nanmean(sd(n_idx).vbap,2), 'filled');
title('VBAP')
xlabel('Azimuth (°)')
ylabel('Elevation (°)')
axis tight
c = colorbar; caxis(lim_colorbar); colormap jet
c.Label.String = 'Spectral distortion (dB)';
set(gca,'FontSize',12)
% set(gca,'Color','k')
filename = [pwd, '\Images\MAP_VBAP_sofaFit2Grid.pdf'];
% exportgraphics(hFigure,filename,'BackgroundColor','none','ContentType','vector')
 
% BILINEAR
hFigure = figure();
scatter(plt_pos(:,1), plt_pos(:,2), 35, nanmean(sd(n_idx).blin,2), 'filled');
title('Bilinear')
xlabel('Azimuth (°)')
ylabel('Elevation (°)')
axis tight
c = colorbar; caxis(lim_colorbar); colormap jet
c.Label.String = 'Spectral distortion (dB)';
set(gca,'FontSize',12)
% set(gca,'Color','k')
filename = [pwd, '\Images\MAP_BILIN_sofaFit2Grid.pdf'];
% exportgraphics(hFigure,filename,'BackgroundColor','none','ContentType','vector')



% Spherical harmonics
hFigure = figure();
scatter(plt_pos(:,1), plt_pos(:,2), 35, nanmean(sd(n_idx).hsph,2), 'filled');
title('Spherical harmonics')
xlabel('Azimuth (°)')
ylabel('Elevation (°)')
axis tight
c = colorbar; caxis(lim_colorbar); colormap jet
c.Label.String = 'Spectral distortion (dB)';
set(gca,'FontSize',12)
% set(gca,'Color','k')
filename = [pwd, '\Images\MAP_HSPH_sofaFit2Grid.pdf'];
% exportgraphics(hFigure,filename,'BackgroundColor','none','ContentType','vector')

%% Visualiazação de posições input e target
n_idx = 4; % 

hFigure = figure();
scatter(plt(n_idx).inpt(:,1), plt(n_idx).inpt(:,2), 28,  [0.3010 0.7450 0.9330], 'filled')
hold on 
scatter(plt(n_idx).des(:,1), plt(n_idx).des(:,2), 28, [0.8500 0.3250 0.0980], 'filled')
title('Input and objective position')
xlabel('Azimuth (°)')
ylabel('Elevation (°)')
axis tight
set(gca,'FontSize',12) 


legend('Input positions',  'Objective position', 'Location', 'best' )
filename = [pwd, '\Images\removedPos_sofaFit2Grid.pdf'];
exportgraphics(hFigure,filename,'BackgroundColor','none','ContentType','vector')


%% Probabilidade distribuição 
% n_idx = 4; % 
% hFigure = figure();
% histogram(mean(sd(n_idx).adpt,2),'Normalization','probability', 'NumBins', 13)
% ylim([0 .35])
% ytix = get(gca, 'YTick');
% set(gca, 'YTick',ytix, 'YTickLabel',ytix*100);
% 
% xlabel('Spectral distortion (dB)')
% ylabel('Probabilidade [%]')
% title('Distribuição da Spectral distortion (Nearest)')
% xlim([0 22.5])
% xticks(0:2:100)
% set(gca,'FontSize',12)
% filename = [pwd, '\Images\Prob_fit2grid_ADPT.pdf'];
% % exportgraphics(hFigure,filename,'BackgroundColor','none','ContentType','vector')
% 
% 
% 
% % VBAP
% hFigure = figure();
% histogram(mean(sd(n_idx).vbap,2),'Normalization','probability', 'NumBins', 13)
% ylim([0 .35])
% ytix = get(gca, 'YTick');
% set(gca, 'YTick',ytix, 'YTickLabel',ytix*100);
% 
% xlabel('Spectral distortion (dB)')
% ylabel('Probabilidade [%]')
% title('Distribuição da Spectral distortion (VBAP)')
% xlim([0 22.5])
% xticks(0:2:100)
% set(gca,'FontSize',12)
% filename = [pwd, '\Images\Prob_fit2grid_VBAP.pdf'];
% % exportgraphics(hFigure,filename,'BackgroundColor','none','ContentType','vector')
% 
% 
% 
% % BILINEAR
% hFigure = figure();
% histogram(mean(sd(n_idx).blin,2),'Normalization','probability', 'NumBins', 13)
% ylim([0 .35])
% ytix = get(gca, 'YTick');
% set(gca, 'YTick',ytix, 'YTickLabel',ytix*100);
% 
% xlabel('Spectral distortion (dB)')
% ylabel('Probabilidade [%]')
% title('Distribuição da Spectral distortion (Bilinear)')
% xlim([0 22.5])
% xticks(0:2:100)
% set(gca,'FontSize',12)
% filename = [pwd, '\Images\Prob_fit2grid_BLIN.pdf'];
% % exportgraphics(hFigure,filename,'BackgroundColor','none','ContentType','vector')
% 
% 
% 
% % Spherical harmonics
% hFigure = figure();
% histogram(mean(sd(n_idx).hsph,2),'Normalization','probability', 'NumBins', 13)
% ylim([0 .35])
% ytix = get(gca, 'YTick');
% set(gca, 'YTick',ytix, 'YTickLabel',ytix*100);
% 
% xlabel('Spectral distortion (dB)')
% ylabel('Probabilidade [%]')
% title('Distribuição da Spectral distortion (Spherical harmonics)')
% xlim([0 22.5])
% xticks(0:2:100)
% set(gca,'FontSize',12)
% filename = [pwd, '\Images\Prob_fit2grid_BLIN.pdf'];
% % exportgraphics(hFigure,filename,'BackgroundColor','none','ContentType','vector')


%% PLOT PER POSITION
n_idx =1; % index for n_desired_pos
fs = 44100;
azi = 240; 
elev = 0;
idx_pos = dsearchn(plt(n_idx).des(:,[1, 2]), [azi, elev]);
N = size(REAL{n_idx}.Data.IR,3);
freq = (0:N/2-1)*fs/N;
ch = 1;

figure()
ff = db(abs(fft(squeeze(REAL{n_idx}.Data.IR(idx_pos, ch, :)), N)));
plot(freq, ff(1:N/2), 'LineWidth', 3.5); hold on
ff = db(abs(fft(squeeze(HSPH{n_idx}.Data.IR(idx_pos, ch, :)), N)));
plot(freq, ff(1:N/2), '--r','LineWidth', 1.4);
ff = db(abs(fft(squeeze(VBAP{n_idx}.Data.IR(idx_pos, ch, :)), N)));
plot(freq, ff(1:N/2), '--','LineWidth', 1.4)
ff = db(abs(fft(squeeze(BLIN{n_idx}.Data.IR(idx_pos, ch, :)), N)));
plot(freq, ff(1:N/2), '--','LineWidth', 1.4)
ff = db(abs(fft(squeeze(ADPT{n_idx}.Data.IR(idx_pos, ch, :)), N)));
plot(freq, ff(1:N/2), '--','LineWidth', 1.4)

legend('Ground truth', 'Spherical harmonics', 'VBAP','Bilinear', 'Nearest', 'Location', 'Best')
arruma_fig('% 4.0f','% 2.0f','ponto')
xlabel('Frequency (Hz)')
ylabel('Amplitude (dB)')
if ch == 1
    title(['Azimuth ' num2str(round(des_pos(idx_pos,1))) '°',...
           '  Elevation ' num2str(round(des_pos(idx_pos,2))) '°', ',  left ear'])
else
    title(['Azimuth ' num2str(round(des_pos(idx_pos,1))) '°',...
           '  Elevation ' num2str(round(des_pos(idx_pos,2))) '°', ',  right ear'])
end
set(gca,'FontSize',12)
axis tight









%% %%%% INTERNAL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Obj, SD, ITD_error, ILD_error]= EvaluateDistortions(Obj_stdy, Obj_real, method)
% methods:  'adapt', 'hybrid', 'vbap', 'bilinear', 'spherical_harmonics'
% des_pos: desired positions to interpolate [azi, ele, rad]
% Obj_stdy: SOFA object with input positions only 
% Obj_real: SOFA object with real postions (to be determined)
% Obj: SOFA object with the interpolated positions
% SD: Spectral distortions between the real and the interpolated HRTFs
% ITD_error: ITD diference between real and interpolated HRTF
% ILD_error: ILD diference between real and interpolated HRTF

%%% Estimar posicoes objetivo a partir das posicoes de entrada
    des_pos = Obj_real.SourcePosition;
    Obj     = sofaFit2Grid(Obj_stdy, des_pos, method); % interpolation
%%% Spectral distortion
    fmin = 500; fmax = 18000;
    SD = sofaSpecDist(Obj, Obj_real, fmin,fmax);
%%% ERRO ITD e ILD 
    [ITD_error, ILD_error] = sofa_ITD_ILD_error(Obj, Obj_real, 'time');
end




function Obj = generate_sofa(positions, hrirs, fs)
    Obj = SOFAgetConventions('SimpleFreeFieldHRIR');
    Obj.Data.IR = hrirs;
    Obj.Data.SamplingRate = fs;
    Obj.SourcePosition = positions;
    Obj = SOFAupdateDimensions(Obj);
end



function plot_mag_map(Obj, Obj_out)
    
    % Check sd 
    fmin = 500; fmax = 18000;
    SD = mean(sofaSpecDist(Obj_out, Obj, fmin,fmax), 'all')
    
    %% PLOTA
    plane1 = 'MagHorizontal';
    plane2 = 'EtcHorizontal';
    
    figure()
    SOFAplotHRTF(Obj,plane1); title(['Reference - ' plane1]);
    axis tight
    xlim([0 2e4])
     
    figure()
    SOFAplotHRTF(Obj_out,plane1); title(['Interpolated - ' plane1]);
    axis tight
    xlim([0 2e4])
    
end


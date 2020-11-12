clear all; clc
% Avaliação de performance entre modelos de interpolação de posições de HRTFs
% Davi R. Carvalho - Abril/2020

%% LOAD 
local = [pwd '\..\Datasets\'];    
path = dir([local 'HUTUBS\pp*_HRIRs_simulated.sofa']);
for k = 1:length(path)
   dataset(k).dados = SOFAload([path(k).folder, '\',path(k).name], 'nochecks');
end


%% General 
% dataset([1:10]) = [];
fs = dataset(1).dados.Data.SamplingRate;
no_posi = size(dataset(1).dados.SourcePosition,1); % numero total de posições


%% 
n=0; cont=0;
no_readpt = 50:100:350;
for k = no_readpt % número de posições "objetivo" (removidas pro teste)
    n=n+1;
    rng(0) % reset random generator
    idx = randperm(no_posi, k); % Index de posições a serem removidas
    for ks = 1:length(dataset)  % numero de individuos sob analise 
        tic
        % HRIRs
        src_IR   = dataset(ks).dados.Data.IR;
        des_IR   = src_IR(idx,:,:); % RI objetivo
        inpt_IR  = src_IR;          
        inpt_IR(idx,:,:) = [];      % RI de entrada

        % Posições  
        inpt_pos = dataset(ks).dados.SourcePosition;
        des_pos  = inpt_pos(idx,:); % Posicoes objetivo
        inpt_pos(idx,:)  = [];      % Posicoes de entrada    
       
        %save para plot
        plt(n).des  = des_pos;
        plt(n).inpt = inpt_pos;

        % SOFA sem posições objetivo (ENTRADA)
        Obj_stdy = SOFAgetConventions('SimpleFreeFieldHRIR');
        Obj_stdy.Data.IR = inpt_IR;
        Obj_stdy.Data.SamplingRate = fs;
        Obj_stdy.SourcePosition = inpt_pos;

        % SOFA com posições objetivo (OBJETIVO)
        REAL = SOFAgetConventions('SimpleFreeFieldHRIR');
        REAL.Data.IR = des_IR;
        REAL.Data.SamplingRate = fs;
        REAL.SourcePosition = des_pos;
        
        % update metadata
        REAL = SOFAupgradeConventions( SOFAupdateDimensions(REAL));    
        Obj_stdy = SOFAupgradeConventions( SOFAupdateDimensions(Obj_stdy));   

        % error
        [ADPT, sd(n).adpt(:,ks), ITD_error(n).adpt(:,ks), ILD_error(n).adpt(:,ks)] = ...
                                        EvaluateDistortions(Obj_stdy, REAL, 'adapt');
        [VBAP, sd(n).vbap(:,ks), ITD_error(n).vbap(:,ks), ILD_error(n).vbap(:,ks)] = ...
                                        EvaluateDistortions(Obj_stdy, REAL, 'vbap');
        [BLIN, sd(n).blin(:,ks), ITD_error(n).blin(:,ks), ILD_error(n).blin(:,ks)] = ...
                                        EvaluateDistortions(Obj_stdy, REAL, 'bilinear');
        [HSPH, sd(n).hsph(:,ks), ITD_error(n).hsph(:,ks), ILD_error(n).hsph(:,ks)] = ...
                                        EvaluateDistortions(Obj_stdy, REAL, 'spherical_harmonics');
        
        clc; disp(['Iteração: ' num2str(n) ' de ' num2str(length(no_readpt)), ....
                  '    Subject: ' num2str(ks), ' de ' num2str(length(dataset))]); 
        cont = cont+1;
        ETA = toc*(length(no_readpt)*length(dataset)-cont)/3600;
        disp(['ETA: ' num2str(ETA),'h' ])      
    end
end

%% SAVE THEM ALL
% save([pwd '\..\DADOS_TREINAMENTO\workspace_Error_sofaFit2Grid.mat'])
clear all
load([pwd '\..\DADOS_TREINAMENTO\workspace_Error_sofaFit2Grid.mat'])


























%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simplicação do erro para comparação (ITD, ILD)
for k = 1:n
    media_itd(k).adpt = mean(mean([ITD_error(k).adpt], 2));
    media_itd(k).vbap = mean(mean([ITD_error(k).vbap], 2));
    media_itd(k).blin = mean(nanmean([ITD_error(k).blin], 2));
%     media_itd(k).hybr = mean(mean([ITD_error(k).hybr], 2));
    media_itd(k).hsph = mean(mean([ITD_error(k).hsph], 2));
    std_itd(k).adpt = std(mean([ITD_error(k).adpt], 2));
    std_itd(k).vbap = std(mean([ITD_error(k).vbap], 2));
    std_itd(k).blin = std(nanmean([ITD_error(k).blin], 2));
%     std_itd(k).hybr = std(mean([ITD_error(k).hybr], 2));
    std_itd(k).hsph = std(mean([ITD_error(k).hsph], 2));
    
    media_ild(k).adpt = mean(mean([ILD_error(k).adpt], 2));
    media_ild(k).vbap = mean(mean([ILD_error(k).vbap], 2));
    media_ild(k).blin = mean(nanmean([ILD_error(k).blin], 2));
%     media_ild(k).hybr = mean(mean([ILD_error(k).hybr], 2));
    media_ild(k).hsph = mean(mean([ILD_error(k).hsph], 2));
    std_ild(k).adpt = std(mean([ILD_error(k).adpt], 2));
    std_ild(k).vbap = std(mean([ILD_error(k).vbap], 2));
    std_ild(k).blin = std(nanmean([ILD_error(k).blin], 2));
%     std_ild(k).hybr = std(mean([ILD_error(k).hybr], 2));
    std_ild(k).hsph = std(mean([ILD_error(k).hsph], 2));

end

figure()
plot([media_itd.adpt]*1e6, 'linewidth', 1.7);hold on
plot([media_itd.vbap]*1e6, 'linewidth', 1.7);
plot([media_itd.blin]*1e6, 'linewidth', 1.7);
% plot([media_itd.hybr]*1e6, 'linewidth', 1.7);
plot([media_itd.hsph]*1e6, 'linewidth', 1.7);
xticks(1:n);
title('ITD error')
xticklabels(no_readpt)
xlabel('Número de posições estimadas')
ylabel('Erro medio (\mus)')
legend('Adaptation', 'VBAP', 'Bilinear',  'Spherical harmonics','location', 'best')


% figure()
% plot([media_ild.adpt], 'linewidth', 1.7);hold on
% plot([media_ild.vbap], 'linewidth', 1.7);
% plot([media_ild.blin], 'linewidth', 1.7);
% % plot([media_ild.hybr], 'linewidth', 1.7);
% plot([media_ild.hsph], 'linewidth', 1.7);
% xticks(1:n);
% title('ild')
% legend('adapt', 'vbap', 'blin',  'hsph','location', 'best')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simplicação do erro para comparação (LSD)
for l = 1:n
    % média do erro de todas as posições para cada individuo
    pos_sd(l).adpt = nanmean(sd(l).adpt);
%     pos_sd(l).hybr = nanmean(sd(l).hybr);
    pos_sd(l).vbap = nanmean(sd(l).vbap);
    pos_sd(l).blin = nanmean(sd(l).blin);
    pos_sd(l).hsph = nanmean(sd(l).hsph);

    % média de todas as posições e todos os individuos
    media_sd(l).adpt = mean(pos_sd(l).adpt);
    stdev(l).adpt    = std(pos_sd(l).adpt); %std entre media de cada indivíduo

%     media_sd(l).hybr = mean(pos_sd(l).hybr);
%     stdev(l).hybr    = std(pos_sd(l).hybr);

    media_sd(l).vbap = mean(pos_sd(l).vbap);
    stdev(l).vbap    = std(pos_sd(l).vbap);

    media_sd(l).blin = mean(pos_sd(l).blin);
    stdev(l).blin    = std(pos_sd(l).blin);
    
    media_sd(l).hsph = mean(pos_sd(l).hsph);
    stdev(l).hsph    = std(pos_sd(l).hsph);
end


%% Media e desvio da SD para todas as posições e todos os individuos 
% padrao entre individuos 
hFigure = figure();
set(0,'DefaultLineLineWidth',1.5)

fig1 = shadedErrorBar(1:n, [media_sd.adpt], [stdev.adpt],{'Color',colors(0),'LineWidth', 1.7},1,0.2,'lin', 0);hold on
fig2 = shadedErrorBar(1:n, [media_sd.vbap], [stdev.vbap],{'Color',colors(1),'LineWidth', 1.7},1,0.2,'lin', 0);
fig3 = shadedErrorBar(1:n, [media_sd.blin], [stdev.blin],{'Color',colors(2),'LineWidth', 1.7},1,0.2,'lin', 0);
% fig4 = shadedErrorBar(1:n, [media_sd.hybr], [stdev.hybr],{'Color',colors(7),'LineWidth', 1.7},1,0.2,'lin', 0);
fig5 = shadedErrorBar(1:n, [media_sd.hsph], [stdev.hsph],{'Color',colors(6),'LineWidth', 1.7},1,0.2,'lin', 0);hold off

%metadata
h = get(gca,'Children');
legendstr = {'', 'Adaptação', '', 'VBAP', '', 'Bilinear', '', 'Spherical Harmonics'};
legend(h([1 3 5 7 ]), legendstr{[ 8 6 4 2]}, 'Location', 'Best')
xticks(1:n)
xticklabels(no_readpt)
title('')
xlabel('Número de posições estimadas')
ylabel('Erro espectral (dB)')
set(gca,'FontSize',12)
ylim([-1 10])
filename = [pwd, '\Images\ShadedError_sofaFit2Grid.pdf'];
% exportgraphics(hFigure,filename,'BackgroundColor','none','ContentType','vector')


%% Mapa de ERRO por posição 
lim_colorbar = [1 8.5];
n_idx = 4; % 
plt_pos = plt(n_idx).des;
plt_inp = plt(n_idx).inpt;

% ADPT
hFigure = figure();
scatter(plt_pos(:,1), plt_pos(:,2), 35, nanmean(sd(n_idx).adpt,2), 'filled'); 
title('Distorção espectral - (Adaptação)')
xlabel('Azimute [°]')
ylabel('Elevação [°]')
axis tight
c = colorbar; caxis(lim_colorbar); colormap jet
c.Label.String = 'Distorção Espectral [dB]';
set(gca,'FontSize',12)
set(gca,'Color','k');
filename = [pwd, '\Images\MAP_ADAPT_sofaFit2Grid.pdf'];
% exportgraphics(hFigure,filename,'BackgroundColor','none','ContentType','vector')
 
%close  VBAP
hFigure = figure();
scatter(plt_pos(:,1), plt_pos(:,2), 35, nanmean(sd(n_idx).vbap,2), 'filled');
title('Distorção espectral logaritmica - (VBAP)')
xlabel('Azimute [°]')
ylabel('Elevação [°]')
axis tight
c = colorbar; caxis(lim_colorbar); colormap jet
c.Label.String = 'Distorção Espectral [dB]';
set(gca,'FontSize',12)
set(gca,'Color','k')
filename = [pwd, '\Images\MAP_VBAP_sofaFit2Grid.pdf'];
% exportgraphics(hFigure,filename,'BackgroundColor','none','ContentType','vector')
 
% BILINEAR
hFigure = figure();
scatter(plt_pos(:,1), plt_pos(:,2), 35, nanmean(sd(n_idx).blin,2), 'filled');
title('Distorção espectral logaritmica - (Bilinear)')
xlabel('Azimute [°]')
ylabel('Elevação [°]')
axis tight
c = colorbar; caxis(lim_colorbar); colormap jet
c.Label.String = 'Distorção Espectral [dB]';
set(gca,'FontSize',12)
set(gca,'Color','k')
filename = [pwd, '\Images\MAP_BILIN_sofaFit2Grid.pdf'];
% exportgraphics(hFigure,filename,'BackgroundColor','none','ContentType','vector')



% SPHERICAL HARMONICS
hFigure = figure();
scatter(plt_pos(:,1), plt_pos(:,2), 35, nanmean(sd(n_idx).hsph,2), 'filled');
title('Distorção espectral logaritmica - (SPHERICAL HARMONICS)')
xlabel('Azimute [°]')
ylabel('Elevação [°]')
axis tight
c = colorbar; caxis(lim_colorbar); colormap jet
c.Label.String = 'Distorção Espectral [dB]';
set(gca,'FontSize',12)
set(gca,'Color','k')
filename = [pwd, '\Images\MAP_HSPH_sofaFit2Grid.pdf'];
% exportgraphics(hFigure,filename,'BackgroundColor','none','ContentType','vector')

%% Visualiazação de posições input e target
hFigure = figure();
scatter(plt(4).inpt(:,1), plt(4).inpt(:,2), 28, 'filled',  'black')
hold on 
scatter(plt(4).des(:,1), plt(4).des(:,2), 28, 'filled', 'red')
title('Posições objetivo dentro do grid original')
xlabel('Azimute [°]')
ylabel('Elevação [°]')
axis tight
set(gca,'FontSize',12) 


legend('Posições de entrada',  'Posições objetivo', 'Location', 'best' )
filename = [pwd, '\Images\removedPos_sofaFit2Grid.pdf'];
% exportgraphics(hFigure,filename,'BackgroundColor','none','ContentType','vector')



%% Probabilidade distribuição 
% n_idx = 4; % 
% hFigure = figure();
% histogram(mean(sd(n_idx).adpt,2),'Normalization','probability', 'NumBins', 13)
% ylim([0 .35])
% ytix = get(gca, 'YTick');
% set(gca, 'YTick',ytix, 'YTickLabel',ytix*100);
% 
% xlabel('Distorção espectral [dB]')
% ylabel('Probabilidade [%]')
% title('Distribuição da distorção espectral (Adaptação)')
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
% xlabel('Distorção espectral [dB]')
% ylabel('Probabilidade [%]')
% title('Distribuição da distorção espectral (VBAP)')
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
% xlabel('Distorção espectral [dB]')
% ylabel('Probabilidade [%]')
% title('Distribuição da distorção espectral (Bilinear)')
% xlim([0 22.5])
% xticks(0:2:100)
% set(gca,'FontSize',12)
% filename = [pwd, '\Images\Prob_fit2grid_BLIN.pdf'];
% % exportgraphics(hFigure,filename,'BackgroundColor','none','ContentType','vector')
% 
% 
% 
% % SPHERICAL HARMONICS
% hFigure = figure();
% histogram(mean(sd(n_idx).hsph,2),'Normalization','probability', 'NumBins', 13)
% ylim([0 .35])
% ytix = get(gca, 'YTick');
% set(gca, 'YTick',ytix, 'YTickLabel',ytix*100);
% 
% xlabel('Distorção espectral [dB]')
% ylabel('Probabilidade [%]')
% title('Distribuição da distorção espectral (Spherical Harmonics)')
% xlim([0 22.5])
% xticks(0:2:100)
% set(gca,'FontSize',12)
% filename = [pwd, '\Images\Prob_fit2grid_BLIN.pdf'];
% % exportgraphics(hFigure,filename,'BackgroundColor','none','ContentType','vector')


%% PLOT PER POSITION
azi = 90; 
elev = 0;
idx_pos = dsearchn(des_pos(:,[1, 2]), [azi, elev]);
N = size(REAL.Data.IR,3);
freq = (0:N/2-1)*fs/N;
ch = 1;

figure()
ff = db(abs(fft(squeeze(REAL.Data.IR(idx_pos, ch, :)), N)));
plot(freq, ff(1:N/2), 'LineWidth', 3.5); hold on
ff = db(abs(fft(squeeze(HSPH.Data.IR(idx_pos, ch, :)), N)));
plot(freq, ff(1:N/2), '--r','LineWidth', 1.4);
ff = db(abs(fft(squeeze(VBAP.Data.IR(idx_pos, ch, :)), N)));
plot(freq, ff(1:N/2), '--','LineWidth', 1.4)
ff = db(abs(fft(squeeze(BLIN.Data.IR(idx_pos, ch, :)), N)));
plot(freq, ff(1:N/2), '--','LineWidth', 1.4)
ff = db(abs(fft(squeeze(ADPT.Data.IR(idx_pos, ch, :)), N)));
plot(freq, ff(1:N/2), '--','LineWidth', 1.4)

legend('Ground truth', 'Spherical harmonics', 'VBAP','Bilinear', 'Nearest', 'Location', 'Best')
arruma_fig('% 4.0f','% 2.0f','virgula')
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
    Obj = sofaFit2Grid(Obj_stdy, des_pos, method); % interpolation
%%% Erro espectral
    fmin = 400; fmax = 20000;
    SD = sofaSpecDist(Obj, Obj_real, fmin,fmax);
%%% ERRO ITD e ILD 
    [ITD_error, ILD_error] = sofa_ITD_ILD_error(Obj, Obj_real, 'time');
end
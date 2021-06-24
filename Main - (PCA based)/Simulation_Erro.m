clear all; clc
% Avaliação de performance do modelo a partir de HRTFs do HUTUBS database
addpath(genpath([pwd, '\..\Functions']));
addpath(genpath([pwd, '\..\DADOS_TREINAMENTO']));

%% LOAD FILES 
% Network
load('DADOS_TREINAMENTO\net_treinada_CIPIC_ARI_ITA_3D3A');
load('DADOS_TREINAMENTO\target_pca_CIPIC_ARI_ITA_3D3A');

%% Load Antropometria 
load('DADOS_TREINAMENTO\input_TUBSIM.mat');
load('DADOS_TREINAMENTO\remove_TUB.mat');

%% Useful  parameters
fs = 44100;
fmin = 250; fmax = 18000;
[no_samples, no_PC, no_directions, no_channels] = size(PCWs);


%% LOAD TEST DATA 
%%% HRIRs ORIGINAIS
local = [pwd, '\..\Datasets\'];
pathtub_sim = dir([local 'HUTUBS\pp*_HRIRs_simulated.sofa']);
[~,idx_tubsim] = natsortfiles({pathtub_sim.name});
pathtub_sim = pathtub_sim(idx_tubsim, :);
pathtub_sim(remove_TUB) = []; %remover individuos sem antropometria

%%% HRIRs GENERICAS
addpath([pwd '\..\Datasets\Generic HRTFs']);
%Fabian
Obj_gen(1).SF = SOFAload('FABIAN_HRIR_measured_HATO_0.sofa', 'nochecks');

% % MIT large pinna
% Obj_gen(2).SF = SOFAload('mit_kemar_large_pinna.sofa');
% % 
% % % MIT normal pinna
% Obj_gen(3).SF = SOFAload('mit_kemar_normal_pinna.sofa');

% Procesamento
for k = 1:length(Obj_gen)
    % Make same grid
    Obj_gen(k).SF = process2unite(Obj_gen(k).S0F, out_pos, fs, fmin, fmax);
end


%% Simulação 
for subj = 1:size(pathtub_sim, 1)% Numero de indivíduos
    clc; disp(['Individuo: ' num2str(subj)])    
    % Simulação de redes 
    for n = 1:no_channels
        for i = 1:no_directions
            X = net_pca{i, n}(anthro(:,subj, n));
            result = mapminmax('reverse',X,PS{i,n});
            DTF_sim(:, i, n) = (coeffs(:, :, i, n) * result) + med_vec2(:, i, n); 
        end
    end    

    %% RECONSTRUÇÃO DE FASE: (fase mínima + ITD)    
    % calcular ITD
    itd_synth_method = 'adapt';
    itd = itd_synthesis(anthro(1,subj,1), anthro(2,subj,1), out_pos, fs, itd_synth_method);
    
    for l = 1:no_directions
        % Aplicando fase minima e de excesso
        offset = 20.4; % valor arbitrario
        [IR_minL, IR_minR] = phase_job(DTF_sim(:, l, 1), DTF_sim(:, l, 2), ...
                                                   itd(l), out_pos(l,:), offset); 
        hrir_final(l, 1, :) = IR_minL;
        hrir_final(l, 2, :) = IR_minR;
    end    
    
    %% Objetos SOFA
    % Simulado 
    Obj_sim = SOFAgetConventions('SimpleFreeFieldHRIR');
    Obj_sim.Data.IR = hrir_final;
    Obj_sim.Data.SamplingRate = fs;
    Obj_sim.SourcePosition = out_pos;
    Obj_sim = SOFAupdateDimensions(Obj_sim);
    
    % Medido Processamento %%% 
    Obj_med = SOFAload([pathtub_sim(subj).folder '\' pathtub_sim(subj).name], 'nochecks');
   % Make same grid
    Obj_med = process2unite(Obj_med, out_pos, fs, fmin, fmax);
    
    %% Erro espectral 
    % Simulado
    sd_sim(:,subj) = sofaSpecDist(Obj_sim, Obj_med, fmin,fmax, 'posi');
    sd_sim_freq(:,subj) = sofaSpecDist(Obj_sim, Obj_med, fmin,fmax, 'freq');
    %%% Erro ITD e ILD 
    [ITDE_sim(:,subj), ILDE_sim(:,subj)] = sofa_ITD_ILD_error(Obj_sim, Obj_med);    

    % HRTFs Genericas
    for k = 1:length(Obj_gen)
        sd_gen(k).gen(:, subj)  = sofaSpecDist(Obj_gen(k).SF, Obj_med, fmin,fmax, 'posi');
        sd_gen(k).freq(:, subj) = sofaSpecDist(Obj_gen(k).SF, Obj_med, fmin,fmax, 'freq');
        %%% Erro ITD e ILD
        [ITDE(k).gen(:,subj), ILDE(k).gen(:,subj)] = sofa_ITD_ILD_error(Obj_gen(k).SF, Obj_med);
    end
end


















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    PLOTS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
save([pwd '\..\DADOS_TREINAMENTO\workspace_Simulation_Erro.mat'])
% clear all; clc
% load([pwd '\..\DADOS_TREINAMENTO\workspace_Simulation_Erro.mat'])

%% Plot erro ESPECTRAL (Todas as posições)
hFigure = figure('Renderer', 'painters', 'Position', [10 10 2000 450]);

% plotar em ordem 
% yitd = cat(1, mean(sd_gen(1).gen),...
%               mean(sd_gen(2).gen),...
%               mean(sd_gen(3).gen),...
%               mean(sd_sim));
% [ysort, idx_sort] = sort(yitd, 1, 'descend');

% N=9;
% C = linspecer(N);
yitd = cat(1, mean(sd_gen(1).gen),...
              mean(sd_sim));
[ysort, ~] = sort(yitd, 1, 'descend');
% hold on 
% for k = 1:size(ysort, 2)
%     for l = 1:size(ysort, 1)
%         if idx_sort(l,k) == 1
%             h(1) = bar(k, ysort(l,k), 'BarWidth', 1, 'FaceAlpha',1, 'FaceColor', 'c');
%         elseif idx_sort(l,k) == 2
%             h(2) = bar(k, ysort(l,k), 'BarWidth', 1, 'FaceAlpha',1, 'FaceColor', C(6,:));
%         elseif idx_sort(l,k) == 3
%             h(3) = bar(k, ysort(l,k), 'BarWidth', 1, 'FaceAlpha',1, 'FaceColor', 'y');
%         elseif idx_sort(l,k) == 4
%             h(4) = bar(k, ysort(l,k), 'BarWidth', 1, 'FaceAlpha',1, 'FaceColor', [0.4940 0.1840 0.5560]);
%         end
%     end
% end
bar(ysort.',  'BarWidth',2 )


% plotar medias 
yline(mean(mean(sd_gen(1).gen)), '--','linewidth', 1.5);
% yline(mean(mean(sd_gen(2).gen)), '--', 'Média', 'color', C(6,:), 'linewidth', 1.5);
% yline(mean(mean(sd_gen(3).gen)), '--y','Média', 'linewidth', 1.5);
yline(mean(mean(sd_sim)),'--', 'Média','linewidth', 1.5); 
hold off

% general
yticks([3.5, mean(mean(sd_sim)), mean(mean(sd_gen(1).gen)),...
          7]);
arruma_fig('% 4.0f','% 2.2f','virgula')
axis([0 95, 3.5, 7.4]);
xlabel('Indivíduo')
ylabel('Erro Espectral [dB]')
title('Distorção espectral por indivíduo')
legend( 'Fabian','Neural network');
set(gca,'fontsize', 12)

filename = [pwd, '\Images\Overall_spectral_error.pdf' ];
% % exportgraphics(hFigure,filename,'BackgroundColor','none','ContentType','vector')


%% Plot ESPECTRAL (plano horizontal) 
hFigure = figure('Renderer', 'painters', 'Position', [10 10 1500 300]);

%Selecionar posições
elev = 0; % Plano horizontal
idx_pos = find(out_pos(:,2)==elev);

% plotar em ordem 
yitd = cat(1, mean(sd_gen(1).gen(idx_pos, :)),...
              mean(sd_gen(2).gen(idx_pos, :)),...
              mean(sd_gen(3).gen(idx_pos, :)),...
              mean(sd_sim(idx_pos, :)));
% yitd = cat(1, mean(sd_gen(1).gen(idx_pos, :)),...
%               mean(sd_sim(idx_pos, :)));
[ysort, idx_sort] = sort(yitd, 1, 'descend');

% Cores
N=9;
C = linspecer(N);

hold on 
for k = 1:size(ysort, 2)
    for l = 1:size(ysort, 1)
        if idx_sort(l,k) == 1
            h(1) = bar(k, ysort(l,k), 'BarWidth', 1, 'FaceAlpha',1, 'FaceColor', 'c');
%         elseif idx_sort(l,k) == 2
%             h(2) = bar(k, ysort(l,k), 'BarWidth', 1, 'FaceAlpha',1, 'FaceColor', C(6,:));
%         elseif idx_sort(l,k) == 3
%             h(3) = bar(k, ysort(l,k), 'BarWidth', 1, 'FaceAlpha',1, 'FaceColor', 'y');
        elseif idx_sort(l,k) == 4
            h(2) = bar(k, ysort(l,k), 'BarWidth', 0.5, 'FaceAlpha',1, 'FaceColor', [0.4940 0.1840 0.5560]);
        end
    end
end

bar(ysort.',  'BarWidth',2 )

% plotar medias 
yline(mean(mean(sd_gen(1).gen(idx_pos, :))), '--', 'Mean', 'linewidth', 1);
yline(mean(mean(sd_gen(2).gen(idx_pos, :))), '--', 'Média', 'color', C(6,:),'linewidth', 1.5);
yline(mean(mean(sd_gen(3).gen(idx_pos, :))), '--', 'color', 'y','linewidth', 1.5);
yline(mean(mean(sd_sim(idx_pos, :))),'--', 'Mean','linewidth', 1); 

hold off


% general 
yticks([3.5, mean(mean(sd_sim(idx_pos, :))), mean(mean(sd_gen(1).gen(idx_pos, :))),...
          7]);
arruma_fig('% 4.0f','% 2.2f','ponto');
axis([0 95, 3.5, 7.4]);
xlabel('Subject');
ylabel('SD (dB)');
legend( 'Fabian','Neural network');
title('Spectral distortion - horizontal plane')
set(gca,'fontsize', 12);

filename = [pwd, '\Images\article_SD_horizontal.pdf' ];
% exportgraphics(hFigure,filename,'BackgroundColor','none','ContentType','vector')



%% MAPA erro espectral por POSICAO
color_range = [1 8];

% Simulada 
figure()
scatter(out_pos(:,1), out_pos(:,2), 30, mean(sd_sim, 2), 'filled', 'square')
title('DTFs Simuladas')
xlabel('Azimute [°]')
ylabel('Elevação [°]')
axis tight
c = colorbar; caxis(color_range); colormap jet
c.Label.String = 'Distorção Espectral [dB]';
set(gca,'FontSize',12)
set(gca,'Color','k')

% export_fig([pwd, '\Images\MAP_sim_simulation'], '-pdf', '-transparent');


% Generica
for i = 1:length(sd_gen) % plot loop
figure()
pl_gen = scatter(out_pos(:,1), out_pos(:,2), 30, mean(sd_gen(i).gen, 2), 'filled', 'square');
switch i 
    case 1
        gen_title = 'DTFs Fabian';
    case 2 
        gen_title = 'DTFs Kemar (orelha maior)';
    case 3 
        gen_title = 'DTFs Kemar (orelha menor)';
end       
title(gen_title)
xlabel('Azimute [°]')
ylabel('Elevação [°]')
axis tight
c = colorbar; caxis(color_range); colormap jet
c.Label.String = 'Distorção Espectral [dB]';
set(gca,'FontSize',12)
set(gca,'Color','k')
% 
% export_fig([pwd, '\Images\MAP_gen_simulation' num2str(i)], '-pdf', '-transparent');
end



%% Histograma da SD por posição 

% hFigure = figure();
% histogram(mean(sd_sim, 2),'Normalization','probability',...
%                                  'NumBins', 17,'orientation','horizontal')
% xlim('auto')
% ylim([2 15])
% xtix = get(gca, 'XTick');
% set(gca, 'XTick',xtix, 'XTickLabel',xtix*100);
% 
% ylabel('Distorção espectral [dB]')
% xlabel('Probabilidade [%]')
% title('Distribuição erro (Simulado)')
% yticks(0:2:100)
% set(gca,'FontSize',12)
% filename = [pwd, '\Images\Prob_SD_sim.pdf'];
% % exportgraphics(hFigure,filename,'BackgroundColor','none','ContentType','vector')
% 
% 
% %%% Genericas
% for i = 1:length(sd_gen) % plot loop
% hFigure = figure();
% histogram(mean(sd_gen(i).gen, 2),'Normalization','probability',...
%                                  'NumBins', 17,'orientation','horizontal')
% xlim('auto')
% xtix = get(gca, 'XTick');
% set(gca, 'XTick',xtix, 'XTickLabel',xtix*100);
% 
% switch i 
%     case 1
%         gen_title = 'Distribuição erro (Fabian)';
%     case 2 
%         gen_title = 'Distribuição erro (Kemar - orelha maior)';
%     case 3 
%         gen_title = 'Distribuição erro (Kemar - orelha menor)';
% end       
% title(gen_title)
% ylabel('Distorção espectral [dB]')
% xlabel('Probabilidade [%]')
% ylim([2 14])
% yticks(0:2:100)
% set(gca,'FontSize',12)
% filename = [pwd, '\Images\Prob_SD_gen_' num2str(i) '.pdf'];
% % exportgraphics(hFigure,filename,'BackgroundColor','none','ContentType','vector')
% end




 %% Plot erro ILD 
% hFigure = figure('Renderer', 'painters', 'Position', [10 10 2000 450]);
% hold on 
% 
% % Selecionar posições
% elev = 0;
% idx_pos = find(out_pos(:,2)==elev);
% 
% 
% % plotar em ordem 
% yild = cat(1, mean(ILDE(1).gen(idx_pos,:)),...
%               mean(ILDE(2).gen(idx_pos,:)),...
%               mean(ILDE(3).gen(idx_pos,:)),...
%               mean(ILDE_sim(idx_pos,:)));
% [ysort, idx_sort] = sort(yild, 1, 'descend');
% 
% % Cores
% N=9;
% C = linspecer(N);
% 
% hold on 
% for k = 1:size(ysort, 2)
%     for l = 1:size(ysort, 1)
%         if idx_sort(l,k) == 1
%             h(1) = bar(k, ysort(l,k), 'BarWidth', 1, 'FaceAlpha',1, 'FaceColor', 'c');
%         elseif idx_sort(l,k) == 2
%             h(2) = bar(k, ysort(l,k), 'BarWidth', 1, 'FaceAlpha',1, 'FaceColor', C(6,:));
%         elseif idx_sort(l,k) == 3
%             h(3) = bar(k, ysort(l,k), 'BarWidth', 1, 'FaceAlpha',1, 'FaceColor', 'y');
%         elseif idx_sort(l,k) == 4
%             h(4) = bar(k, ysort(l,k), 'BarWidth', 1, 'FaceAlpha',1, 'FaceColor', [0.4940 0.1840 0.5560]);
%         end
%     end
% end
% 
% % plotar medias 
% yline(mean(mean(ILDE(1).gen(idx_pos,:))), '--', 'Média', 'color', 'c', 'linewidth', 1.5);
% yline(mean(mean(ILDE(2).gen(idx_pos,:))), '--', 'Média', 'color', C(6,:),'linewidth', 1.5);
% yline(mean(mean(ILDE(3).gen(idx_pos,:))), '--', 'color', 'y','linewidth', 1.5);
% yline(mean(mean(ILDE_sim(idx_pos,:))),'--', 'Média', 'color', [0.4940 0.1840 0.5560],'linewidth', 1.5); 
% hold off
% 
% % general 
% yticks(sort([3, mean(mean(ILDE_sim(idx_pos,:))), mean(mean(ILDE(3).gen(idx_pos,:))),...
%            mean(mean(ILDE(1).gen(idx_pos,:))), mean(mean(ILDE(2).gen(idx_pos,:))), 8]));
% arruma_fig('% 4.0f','% 2.2f','virgula')
% axis tight
% title('Erro médio do ILD no plano horizontal')
% xlabel('Indivíduo')
% ylabel('Erro Médio Absoluto')
% legend(h, 'Fabian','Kemar Small','Kemar Large','Simulada','location','southwest')
% set(gca,'fontsize', 12)
% 
% filename = [pwd, '\Images\Overall_ILD_error.pdf' ];
% % % exportgraphics(hFigure,filename,'BackgroundColor','none','ContentType','vector')
% 

%% Plot polar ITD
% Selecionar posições
elev = 0;
idx_pos = find(round(out_pos(:,2))==elev);

figure('Renderer', 'painters', 'Position', [10 10 700 450])

%Cores
N=5;
C = linspecer(N);

angles = deg2rad(out_pos(idx_pos, 1));
ax = polaraxes;
polarplot(angles, smoothdata(mean(ITDE(2).gen(idx_pos,:),2)*10^6), ...
                  '.-', 'color', 'k', 'linewidth', 1.2); hold on
polarplot(angles, smoothdata(mean(ITDE(3).gen(idx_pos,:),2)*10^6), ...
                  '.-',  'color', C(3,:) , 'linewidth', 1.2);
polarplot(angles, smoothdata(mean(ITDE(1).gen(idx_pos,:),2)*10^6),...
                  '.-','color', C(1,:), 'linewidth', 1.2); hold on
polarplot(angles, smoothdata(mean(ITDE_sim(idx_pos,:),2)*10^6), ...
                  '.-','color', C(2,:), 'linewidth', 1.2); hold off 
ax.ThetaDir = 'counterclockwise';
ax.ThetaZeroLocation = 'top';
% title('Erro médio do ITD')
legend('Kemar Large','Kemar Small', 'Fabian', 'Simulado: Esférico')
% legend('Fabian', 'Neural network')

thetaticks(0:30:330)
thetaticklabels({'0°', '30°', '60°', '90°', '120°', '150°', '180°', '210°', '240°',...
                 '270°', '300°', '330°'})
rticks([50 100 150])
rticklabels({'50 \mus','100 \mus','150 \mus'})

axis tight
grid on
ax.GridAlpha = 0.3; 
set(gca,'fontsize', 12)
% % export_fig([pwd, '\Images\Erro_ITD_horizontal_SPHERIC' ], '-pdf', '-transparent');

%% Plot erro ITD 
% % figure('Renderer', 'painters', 'Position', [10 10 2000 450])
% % bar(mean(ITDE(1).gen),'BarWidth', 1, 'FaceAlpha',0.5, 'FaceColor','c'); hold on
% % bar(mean(ITDE(2).gen),'BarWidth', 1, 'FaceAlpha',0.5, 'FaceColor','k');
% % bar(mean(ITDE(3).gen),'BarWidth', 1, 'FaceAlpha',0.5, 'FaceColor','y'); 
% % bar(mean(ITDE_sim),'BarWidth', 1, 'FaceAlpha',0.5, 'FaceColor',[0.4940 0.1840 0.5560]); 
% % 
% % 
% % % plot again to make the smallest in front
% % yitd = cat(1, mean(ITDE(1).gen), mean(ITDE(2).gen), mean(ITDE(3).gen), mean(ITDE_sim));
% % [~, rol] = min(yitd);
% % for k = 1:length(yitd)
% %     if rol(k) == 1
% %             yit = mean(ITDE(1).gen);
% %             bar(k, yit(k), 'BarWidth', 1, 'FaceAlpha',0.6, 'FaceColor','c')
% %     elseif rol(k) == 2
% %             yit = mean(ITDE(2).gen);
% %             bar(k, yit(k), 'BarWidth', 1, 'FaceAlpha',0.6, 'FaceColor','k')   
% %     elseif rol(k) == 3
% %             yit = mean(ITDE(3).gen);
% %             bar(k, yit(k), 'BarWidth', 1, 'FaceAlpha',0.6, 'FaceColor','y')
% %     elseif rol(k) == 4
% %             yit = mean(ITDE_sim);
% %             bar(k, yit(k), 'BarWidth', 1, 'FaceAlpha',0.6, 'FaceColor',[0.4940 0.1840 0.5560])
% %     end
% % end
% % % axis([1 93, 2.3, 7.4])
% % legend(  'Fabian', 'Kemar large', 'Kemar small','Simulada', 'location', 'southeast')
% % title('ITD')
% % xlabel('Indivíduo')
% % ylabel('Erro Médio Absoluto')
% % set(gca,'fontsize', 12)
% % arruma_fig('% 4.0f','% 2.2f','virgula')
% 
% % % export_fig([pwd, '\Images\Overall_ITD_error' ], '-pdf', '-transparent');
% 
% 

%% LOCAL FUCTIONS 
function Obj = process2unite(Obj, out_pos, fs, fmin, fmax)
%     % Make same grid
    Obj = sofaFit2Grid(Obj, out_pos, 'Fs', fs);     
    % filter
    Obj = sofaIRfilter(Obj, fmin, fmax);
    % HRTF -> DTF
    [Obj, ~] = SOFAhrtf2dtf(Obj);
    % Normalize L/R balance and IR levels
    Obj = sofaNormalize(Obj);
end
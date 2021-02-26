% clear all; close all; clc
% DAVI ROCHA CARVALHO @ENG. ACUSTICA - UFSM; MARÇO/2020 
% Avaliacao da IR simulada, por meio de comparacoes entre IR medida,
% simulada e generica
% Comparacao do espectro, ILD e ITD
addpath([pwd, '\..\Functions']); 

%% Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defina quais datasets usar: {'cipic', 'ari', 'ita', '3d3a', 'riec', 'tub'}

Datasets = {'cipic', 'ari', 'ita', '3d3a'};

% Range de frequencias de ANALISE
fmin = 250; fmax = 18000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Simulated and Measured SOFAs
file = '\..\DADOS_TREINAMENTO\rebuilt_data';
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
load([pwd file])
disp('Dados carregados!')

%% 
fs = Obj_med.Data.SamplingRate;
out_pos = Obj_med(1).SourcePosition;

%% LOAD Generic HRTFs]
addpath([pwd, '\..\Datasets\Generic HRTFs\']);
%Fabian
generic_head(1).SF = SOFAload('FABIAN_HRIR_measured_HATO_0.sofa');

% MIT large pinna
generic_head(2).SF = SOFAload('mit_kemar_large_pinna.sofa');

% MIT normal pinna
generic_head(3).SF = SOFAload('mit_kemar_normal_pinna.sofa');

% Aplicar mesmo pre processamento feito na HRIR simulada e medida
for k = 1:length(generic_head)
    Obj_gen(k).SF = sofaFit2Grid(generic_head(k).SF, out_pos,'adapt', 'Fs', fs);
    Obj_gen(k).SF = sofaIRfilter(Obj_gen(k).SF, fmin, fmax); 
    Obj_gen(k).SF = SOFAhrtf2dtf(Obj_gen(k).SF);
end

%%
% export_fig([pwd, '\Images\pos_distance_kemar'], '-pdf', '-transparent');

%% Logarithmic Spectral Distortion (LSD)
sd_sim = sofaSpecDist(Obj_sim, Obj_med, fmin,fmax);
for l = 1:length(generic_head)
    sd_gen(:,l) = sofaSpecDist(Obj_gen(l).SF, Obj_med, fmin,fmax);
end

%% Mapa projeção
% color_range = [1 8.5];
% 
% % Simulada 
% figure()
% scatter(out_pos(:,1), out_pos(:,2), 30, sd_sim, 'filled', 'square')
% title('HRTFs Simuladas')
% xlabel('Azimute [grau]')
% ylabel('Elevação [grau]')
% axis tight
% c = colorbar; caxis(color_range); colormap jet
% c.Label.String = 'Distorção Espectral [dB]';
% set(gca,'FontSize',12)
% set(gca,'Color','k')
% 
% % export_fig([pwd, '\Images\MAP_hrtf_sim'], '-pdf', '-transparent');
% 
% 
% % Generica
% for i = 1:length(generic_head) % plot loop
% figure()
% pl_gen = scatter(out_pos(:,1), out_pos(:,2), 30, sd_gen(:, i), 'filled', 'square');
% switch i 
%     case 1
%         gen_title = 'HRTFs Fabian';
%     case 2 
%         gen_title = 'HRTFs Kemar (orelha maior)';
%     case 3 
%         gen_title = 'HRTFs Kemar (orelha menor)';
% end       
% title(gen_title)
% xlabel('Azimute [grau]')
% ylabel('Elevação [grau]')
% axis tight
% c = colorbar; caxis(color_range); colormap jet
% c.Label.String = 'Distorção Espectral [dB]';
% set(gca,'FontSize',12)
% set(gca,'Color','k')
% 
% % export_fig([pwd, '\Images\MAP_hrtf_gen_' num2str(i)], '-pdf', '-transparent');
% end


%% ITD & ILD error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:length(generic_head)
    [ITDE_gen(:,k), ILDE_gen(:,k)] = sofa_ITD_ILD_error(Obj_gen(k).SF, Obj_med);
end

[ITDE_sim, ILDE_sim] = sofa_ITD_ILD_error(Obj_sim, Obj_med);


%% Plot polar ITD
% Erro ITD no plano horizontal
elev = 0;
idx_pos = find(out_pos(:,2)==elev);

angles = deg2rad(out_pos(idx_pos, 1));

figure('Renderer', 'painters', 'Position', [10 10 700 450])
N=5;
C = linspecer(N);
ax = polaraxes;
polarplot(angles, smoothdata(mean(ITDE_gen(idx_pos,2),2)*10^6), ...
                  '.-', 'color', 'k', 'linewidth', 1.2); hold on
polarplot(angles, smoothdata(mean(ITDE_gen(idx_pos,3),2)*10^6), ...
                  '.-',  'color', C(3,:) , 'linewidth', 1.2);
polarplot(angles, smoothdata(mean(ITDE_gen(idx_pos,1),2)*10^6),...
                  '.-','color', C(1,:), 'linewidth', 1.2); hold on
polarplot(angles, smoothdata(mean(ITDE_sim(idx_pos,:),2)*10^6), ...
                  '.-','color', C(2,:), 'linewidth', 1.2); hold off 

ax.ThetaDir = 'counterclockwise';
ax.ThetaZeroLocation = 'top';
% title('Erro médio do ITD')
legend( 'Kemar large','Kemar small',  'Fabian', 'Simulated')

thetaticks(0:30:330)
thetaticklabels({'0°', '30°', '60°', '90°', '120°', '150°', '180°', '210°', '240°',...
                 '270°', '300°', '330°'})
rticks([50 100 150])
rticklabels({'50 \mus','100 \mus','150 \mus'})

axis tight
grid on
ax.GridAlpha = 0.3; 
set(gca,'fontsize', 12)



%% GENERAL ERROR
clc
% Espectro
media_LSD_gen = mean(sd_gen)
media_LSD_sim = mean(sd_sim)
% Tempo 
ITDE_gen_mu = mean(ITDE_gen)*10^6
ITDE_sim_mu = mean(ITDE_sim)*10^6
% Amplitude
ILDE_gen_mu = rms(ILDE_gen)
ILDE_sim_mu = rms(ILDE_sim)

% Media geral
media_erro_sim = (media_LSD_sim );
media_erro_gen = (media_LSD_gen );

% figure()
% hold on 
% bar(media_erro_sim,'g','FaceAlpha',0.8); 
% bar(2:4, media_erro_gen,'r','FaceAlpha',0.8);
% hold off
% legend('HRIRs Simuladas', 'HRIRs Genericas', 'Location', 'northwest');
% title('Média do erro espectral, itd e ild')
% xticks([1:4]);
% xticklabels({'Simulada', 'Fabian', 'Kemar large', 'Kemar small'})
% ylabel('Erro espectral [dB]')

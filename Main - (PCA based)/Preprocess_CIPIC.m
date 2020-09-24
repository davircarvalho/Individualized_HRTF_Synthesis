clear all; clc; 
% DAVI ROCHA CARVALHO; ENG. ACUSTICA - UFSM; MAIO/2019
% Analise antropométrica inicial, definição de parâmetros input para ANN
% por meio de Principal Component Analysis (PCA) - CIPIC dataset
addpath(genpath([pwd, '\..\EAC-Toolbox']));

tic
%% Carregar o banco de HRIRs 
% Defina como path o local onde estão as pastas com os arquivos .mat do
% banco de dados.
onde = 'B:\Documentos\#3 - TCC\EAC-TCC-Davi\HRTF-datasets\Banco CIPIC\standard_hrir_database';
path = dir([onde '\subject_*\hrir_final.mat']);
% Importa todo o banco de indivíduos para 'data'
for jj = 1 : length( path )
    data(jj) = importdata( fullfile( path(jj).folder, path(jj).name) );
end

% Remove indivíduos sem dados antropométricos e 1 indivíduo pra teste MANUAL
% individuos teste: 1 (id:'003')
remove = [1:3, 5:8, 10, 42];
fields = fieldnames(data); %salva nome dos fields do struct para reconstrução 
data  = struct2cell(data); %transforma em celula para remover indivíduos
data(:,:, remove) = [];    %remove indivíduos 
data = cell2struct(data, fields, 1); %volta para struct



%% PLOT same position multiple subjects
% % widescreen
% % figure('Renderer', 'painters', 'Position', [10 10 2000 562])
% 
% for k = 1 : length(data)
%     
%    y(:, k) = db(abs(fft(squeeze(data(k).hrir_l(13, 9, :)))));
%     
% end
% N = 200;
% fs = 44100;
% freq = linspace(0, fs-fs/N, N);
% plot(freq(1:N/2), y(1:N/2, :), 'linewidth', 1.2)
% 
% % % PT_BR
% xlabel('Frequência [Hz]')
% ylabel('Amplitude [dB]')
% title('HRTFs orelha esquerda (azimute 0°, elevação 0°)')
% % labels = {'Individuo 03','Individuo 08','Individuo 09','Individuo 10', 'Individuo 11'};
% 
% % % EN_US
% % xlabel('Frequency [Hz]')
% % ylabel('Amplitude [dB]')
% % title('HRTFs left ear (azimuth 0°, elevation 0°)')
% % labels = {'Subject 03','Subject 08','Subject 09','Subject 10', 'Subject 11'};
% 
% 
% % legend(labels, 'location', 'best')
% set(gca, 'FontSize', 13);
% arruma_fig('% 4.0f','% 2.0f','virgula')
% xlim([800 20000])
% ylim([-30, 14])
% grid on
% % set(gca, 'YLimSpec', 'Tight');
% % export_fig([pwd, '\Images\English\intro_hrtfs' ], '-pdf', '-transparent');
% export_fig([pwd, '\Images\intro_hrtfs' ], '-pdf', '-transparent');



%% Extração de dados [HRIR]
% Veja o grid de dispsição das caixas
no_subjects = length(data);
for jj = 1:no_subjects
    hrir(jj).L = data(jj).hrir_l; 
    hrir(jj).R = data(jj).hrir_r;
    
    temp = data(jj).ITD; 
    %transforma a matriz em um vetor coluna, em que todos os azimutes 
    % para 1 elevação são postos em sequencia.
    ITD_cipic(:, jj) = reshape(temp,[], 1); 
end


%% HRTF EM LOG SCALE [HRTF_log]
% Aplica fft e transforma para escala log
no_samples = length(hrir(1).L);
for jj = 1: no_subjects
    for m = 1:25
        for n = 1:50
            hrtf(jj).logL(m,n,:) = 20*log10(abs(fft(hrir(jj).L(m,n,:), no_samples))./ no_samples); 
            hrtf(jj).logR(m,n,:) = 20*log10(abs(fft(hrir(jj).R(m,n,:), no_samples))./ no_samples);
        end
    end
end
[dim1, dim2, no_samples] = size(hrtf(1).logL);

%% [HRTF_dir_log]

% Fazer a média de TODAS AS HRTFS DE UM INDIVIDUO e em seguida subtraimos a
% média de cada hrtf a fim de remover componentes que independem da
% direção
% A média é conhecida também como CTF (comum transfer function)
% A magnitude subtraida da média é conhecida como DTF (directional transfer functions)
no_directions = dim1*dim2; %1250 direções
no_channels = 2;

% somatório de todas as direções para determinado indivíduo e orelha
% calculating directional mean
med_vec = zeros(no_samples, no_subjects, no_channels);

for i = 1:no_subjects
    temp1 = zeros(no_samples,1); 
    for j = 1:dim1
       for k = 1:dim2
            temp1 = temp1 + squeeze(hrtf(i).logL(j, k, :)) + squeeze(hrtf(i).logR(j, k, :));
       end
    end
    temp1 = temp1/(2*no_directions);   
    med_vec(:, i) = temp1; % (no_samples x no_subjecs x no_channels)
end

% "the mean is subtracted from each HRTFlog to obtain a new transfer
% function (HRTF_dir_log)
%  which represents primarily direction-dependent spectral effects"

for i = 1: no_subjects
    for j = 1:dim1
        for k = 1:dim2
            hrtf(i).dir_logL(j, k, :) = squeeze(hrtf(i).logL(j, k, :)) - med_vec(:, i);
            hrtf(i).dir_logR(j, k, :) = squeeze(hrtf(i).logR(j, k, :)) - med_vec(:, i);
        end
    end
end


%% Reestruturação dos dados para input na PCA
for jj = 1:no_subjects
    pos = 1;
    for k = 1:dim1
        for l = 1:dim2         
            HRTF_d(:, pos, jj, 1) = hrtf(jj).dir_logL(k, l, 1:no_samples/2);
            HRTF_d(:, pos, jj, 2) = hrtf(jj).dir_logR(k, l, 1:no_samples/2);
            pos = pos + 1;
        end
    end
end
DTF_CIPIC = permute(HRTF_d, [1, 3, 2, 4]);

%% Smooth da PTF
% Josef recomenda aplicar smooth nos dados antes de aplicar a PCA, melhora
% a compressão nos dados por reduzir a variabilidade nos dados sem diminuir
% a precisão na localização
for k = 1:no_channels
    for l = 1:no_directions
        for m = 1:no_subjects
            DTF_CIPIC(:,m,l,k) = smooth(DTF_CIPIC(:,m,l,k), 'sgolay', 3);
        end
    end
end

 %% Principal Component Analysis (PCA) - METODO 1
no_PC = 15; %número de principais componentes de interesse
% [weight_vectors, PC_mtx2, eigen_value] = PCA_fun(PTF, no_PC);
% 12x 30x 1250x 2 - formato de input na neural net 

% Para calcular o quanto representa o numero de princcipais componentes
% representa o dataset inteiro, assumimos que com todos as PCs teremos 100%
% logo para normalizar é necessário somar todos os valores  de variancia e dividir o
% vetor original pela soma, assim ao somarmos os n primeiros valores,
% teremos a porcentagem do quanto aquela quantidade de PC representa o dataset original

% A PCA é aplicada para todos os sujeitos em uma mesma direção
for m = 1:no_channels
    for n = 1:no_directions        
        med_vec2(:,n,m) = mean(DTF_CIPIC(:, :, n, m),2);
        data_mtx = DTF_CIPIC(:, :, n, m) - med_vec2(:,n,m); % centralizar PCA
         
        % Principal component analysis
        [coeff, score, latent] = pca(data_mtx, 'NumComponents', no_PC, ...
                                                  'Centered', false);
        %  "Component scores are the transformed variables
        % and loadings(coeff) describe the weights to multiply the normed original variable
        % to get the component score."
        % In this case pca does not center the data.
        % You can reconstruct the original data using score*coeff'              
        coeffs(:,:, n, m) = coeff; %Cada matrix coeff tem que ser orthonormal
        PCWs(:,:, n, m) = score;
        variancia(:, n, m) = latent;
    end
end
disp('PCA calculada!')


%% PLOT PCA RECONSTRUCTION PER PC 
close all
dir = 150; % direção 
npc = 9;   %  numero de cp
subj = 10; % sujeito
ch = 1;    % orelha

recon = PCWs(:,1:npc, dir, ch)*coeffs(:,1:npc,dir,ch)' + med_vec2(:,dir,ch); % for all subj

% vetor de frequencia 
N = 200; fs = 44100;
freq = (0:N-1)/N*fs;
fr = freq(1:N/2);

figure()
semilogx(fr, DTF_CIPIC(:,subj,dir,ch), 'linewidth', 2); hold on
semilogx(fr, recon(:,subj), '--r', 'linewidth', 1.5); hold off
legend({'Original', ['Reconstrução: ' num2str(npc) ' CP']}, 'location', 'best')
axis tight
xlabel('Frequência [Hz]')
ylabel('Amplitude [dB]')
arruma_fig('% 4.0f','% 2.0f','virgula')
set(gca, 'FontSize', 13);
export_fig([pwd, '\Images\pc_recon_' num2str(npc) '_CP'], '-pdf', '-transparent');
%% Dados ANTHROPOMETRICOS
% Carregar dados 
anthro1 = load('anthro_CIPIC.mat');

% separando os 8 parametros necessários para a ann
% d1, d3 : d6
d1 = anthro1.D(:, [1, 3:6]);   % L
d2 = anthro1.D(:, ([9, 11:14])); % R
% x1, x3, x12
x = anthro1.X(:, [1, 3, 12]);

% unir as duas matrizes
anthro(:,:,1) = [d1, x]';  % L
anthro(:,:,2) = [d2, x]';  % R

%removeremos os individuos que não foram medidos e indivíduos para teste
% indivíduos para teste: 1, 4, 9, 11:14.
anthro(:,remove,:) = [];

%% SALVAR DADOS
save('DADOS_TREINAMENTO\input_cipic.mat','anthro');
save('DADOS_TREINAMENTO\target_pca_cipic.mat', 'PCWs', 'coeffs', 'med_vec2');

disp('Dados salvos!')
toc

    
%% PLOT DO NUMERO DE PC NECESSÀRIOS
% Para calcular o quanto representa o numero de princcipais componentes
% representa o dataset inteiro, assumimos que com todos as PCs teremos 100%
% logo para normalizar é necessário somar todos os valores  de variancia e dividir o
% vetor original pela soma, assim ao somarmos os n primeiros valores,
% teremos a porcentagem do quanto aquela quantidade de PC representa o dataset original

% numero de PC vs variancia
vari = variancia(:,1,1);
pn_pc = 1:length(vari(:,1,1));
var_norm = vari./sum(vari);
for i = 1:length(vari)
    var_var(i) = sum(var_norm(1:i));
end

figure()
plot(pn_pc, var_var*100, 'LineWidth', 2.0)
axis tight
ylabel('Recuperação de variancia no dataset [%]')
xlabel('Número de Componentes Principais [~]')
grid on 
arruma_fig('% 4.0f','% 2.0f','virgula')
set(gca, 'FontSize', 13);
% xtick = var_var*100;
% xticks 
export_fig([pwd, '\Images\no_PC_cipic' ], '-pdf', '-transparent');



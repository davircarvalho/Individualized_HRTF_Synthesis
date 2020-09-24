clear all; close all; clc;
% DAVI ROCHA CARVALHO; ENG. ACUSTICA - UFSM; FEVEREIRO/2020 
addpath(genpath([pwd, '\..\EAC-Toolbox']));

%% Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defina quais datasets usar: {'cipic', 'ari', 'ita', '3d3a', 'riec', 'tub'}

Datasets = {'cipic', 'ari', 'ita', '3d3a'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD
target_file = 'DADOS_TREINAMENTO\target_pca';
input_file  = 'DADOS_TREINAMENTO\input';

if any(strcmp('cipic', Datasets)) 
    target_file = append(target_file, '_CIPIC');
    input_file  = append(input_file, '_CIPIC');
end
if any(strcmp({'ari'}, Datasets)) 
    target_file = append(target_file, '_ARI');
    input_file  = append(input_file, '_ARI');
end
if any(strcmp({'ita'}, Datasets)) 
    target_file = append(target_file, '_ITA');
    input_file  = append(input_file, '_ITA');
end
if any(strcmp({'3d3a'}, Datasets)) 
    target_file = append(target_file, '_3D3A');
    input_file  = append(input_file, '_3D3A');
end
if any(strcmp({'riec'}, Datasets)) 
    target_file = append(target_file, '_RIEC');
    input_file  = append(input_file, '_RIEC');
end
if any(strcmp({'tub_meas'}, Datasets)) 
    target_file = append(target_file, '_TUBMEAS');
    input_file  = append(input_file, '_TUBMEAS');
end
if any(strcmp({'tub_sim'}, Datasets)) 
    target_file = append(target_file, '_TUBSIM');
    input_file  = append(input_file, '_TUBSIM');
end
load(target_file)
load(input_file)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
target = coeffs;
[no_PC, no_subjects, no_directions, no_channels] = size(target);

%% Pré-processamento INPUT (normalize)
% [xinpt, sig1, mu1] = nn_preprocess(anthro);

%% Setting up Feed Forward Neural Network with BP
net = fitnet(14, 'trainbr'); 
% função de treinamento
% max iterações 
net.trainParam.epochs = 400; 
% max validation non improvement
net.trainParam.max_fail = 4; % avoid overfitting by early stopping

% dividir minibatchs
net.divideFcn= 'divideint'; 
net.divideParam.trainRatio = .8; 
net.divideParam.valRatio   = .15; 
net.divideParam.testRatio  = .05; 

net.layers{1}.transferFcn  = 'tansig';
net.layers{2}.transferFcn  = 'purelin';

net.trainParam.showWindow  = 0; %show training window
%% Treinamento 
disp('Iniciado o treinamento da rede neural')
tic
wait = waitbar(0,'Progresso de treinamento');
cont = 0;
for channel = 1:no_channels % cada orelha
    x = ((anthro(:, :, channel)));  
    for i = 1:no_directions
        cont = cont+1;
        waitbar(cont/(no_directions*no_channels), wait);        
        t = target(:, :, i, channel)';% output target varia com a direção 
        net = configure(net, x, t); %define input e output da rede
               %% VIEW net export
%                 jframe = view(net);
%                 %# create it in a MATLAB figure
%                 hFig = figure('Menubar','none', 'Position',[100 100 565 166]);
%                 jpanel = get(jframe,'ContentPane');
%                 [~,h] = javacomponent(jpanel);
%                 set(h, 'units','normalized', 'position',[0 0 1 1])
%                 %# close java window
%                 jframe.setVisible(false);
%                 jframe.dispose();           
%                 export_fig([pwd, '\Images\view_net' ], '-pdf', '-transparent');
        %% 
            [net_pca{i, channel}, tr_pca{i, channel}] = train(net, x, t);
  
        
%%%%%%% Repete o treinamento para mesma direção em caso de perf muito ruim
        cont1 = 0;
        tol_max = 0.01;
         while tr_pca{i, channel}.best_tperf > tol_max
            net = init(net);           
            [net_pca{i, channel}, tr_pca{i, channel}] = train(net, x, t);
            cont1 = cont1 + 1;
            if cont1 > 10 % caso não consiga atingir a meta em 5 tentativas, utilizará o ultimo treinamento realizado 
                tr_pca{i, channel}.best_tperf
                break
            end                   
         end  
    end
end
disp('Treinamento completo.')
toc
close(wait)
%% SAVE DATA
path_save = 'DADOS_TREINAMENTO\net_treinada';
if any(strcmp('cipic', Datasets))
    path_save = append(path_save, '_CIPIC');
end
if any(strcmp('ari', Datasets))
    path_save = append(path_save, '_ARI');
end
if any(strcmp('ita', Datasets))
    path_save = append(path_save, '_ITA');
end
if any(strcmp('3d3a', Datasets))
    path_save = append(path_save, '_3D3A');
end
if any(strcmp('riec', Datasets))
    path_save = append(path_save, '_RIEC');
end
if any(strcmp('tub_sim', Datasets))
    path_save = append(path_save, '_TUBSIM');
end
if any(strcmp('tub_meas', Datasets))
    path_save = append(path_save, '_TUBMEAS');
end
save(path_save, 'net_pca', 'tr_pca');

disp('Dados salvos!')










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT - MSE for the test set
figure('Renderer', 'painters', 'Position', [10 10 2000 300])
for k = 1:no_channels
    for l = 1:no_directions
        pv(l, k) = tr_pca{l, k}.best_vperf;
    end
end
   
plot(pv(:,1)); hold on
plot(pv(:,2)); hold off
xlim([0 no_directions])
ylim([4e-3 10e-3])

% legend('Esquerda', 'Direita');
% xlabel('Índice da posição')
legend('Left', 'Right');
xlabel('Neural Nets')
ylabel('RMSE')
% title('Performance no grupo de validação')
title('Performance on validation dataset')

set(gca, 'FontSize', 16)

% export_fig([pwd, '\Images\English\neuralnet_vRMSE' ], '-pdf', '-transparent');


%% PLOT - MSE for the training set

figure('Renderer', 'painters', 'Position', [10 10 2000 300])
for k = 1:no_channels
    for l = 1:no_directions
        pt(l, k) = tr_pca{l, k}.best_perf;
    end
end
   
plot(pt(:,1)); hold on
plot(pt(:,2)); hold off
xlim([0 no_directions])
ylim([4e-3 10e-3])
% legend('Esquerda', 'Direita');
% xlabel('Índice da posição')
legend('Left', 'Right');
xlabel('Neural Nets')
ylabel('RMSE')
title('Performance on training dataset')
set(gca, 'FontSize', 16)
% export_fig([pwd, '\Images\English\neuralnet_tRMSE' ], '-pdf', '-transparent');

%% 
% tperf=[];
% for l = 1:size(tr_pca, 1)
%     ttemp = tr_pca{l, 1}.tperf;
%     if length(ttemp) > length(tperf)
%         tperf = [tperf zeros(1,abs(length(tperf) - length(ttemp)))];
%     else %tperf >ttemp
%         ttemp = [ttemp zeros(1,abs(length(tperf) - length(ttemp)))];
%     end
%     tperf = tperf +ttemp;
% end
% 

clear all; close all; clc;
% DAVI ROCHA CARVALHO; ENG. ACUSTICA - UFSM; MAIO/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD
addpath([pwd '/DADOS_TREINAMENTO'])
load('target_pca_cipic.mat')
load('input_cipic')

target = coeffs;
[no_PC, no_subjects, no_directions, no_channels] = size(target);

%% Pré-processamento INPUT TARGET (normalize)
% [xinpt, sig, mu] = nn_input_preprocess(InputMatrix);

%% parameters count 
% [I, ~] = size(InputMatrix_CIPIC);
% [O, Ntrn] = size(weights);
% 
% H = 2*I+1;                   % No of hidden layers
% Nw = (I+1)*H+(H*1)*O;    % No of unknown weights
% Ntrneq = Ntrn*O;         % No of training equations
% Ndof   = Ntrneq - Nw;    % No of degrees of freedom

%%
net = fitnet(11, 'traincgb'); 
net.trainParam.max_fail = 6;
net.trainParam.showWindow = 0; %show training window
net.performFcn='msereg';

net.divideFcn= 'dividerand'; % divide the data randomly 
net.divideParam.trainRatio= .75; % we use 70% of the data for training 
net.divideParam.valRatio= .2; % 30% is for validation
net.divideParam.testRatio= .05; % 0% for testing

tol_max = 0.04; % tolerancia máxima do erro do dataset de testes

%% Treinamento 
disp('Iniciado o treinamento da rede neural')
tic
for channel = 1:no_channels % cada orelha
    x = (anthro(:, :, channel));  
    clc;
    for i = 1:no_directions
        disp(i)
        t = (target(:, :, i, channel));% output target varia com a direção 
        
        net = configure(net, x, t'); %define input e output da rede
        [net_pca{i, channel}, tr_pca{i, channel}] = train(net, x, t');
  
%       Repete o treinamento para mesma direção em caso de perf muito ruim
        cont1 = 0;
         while tr_pca{i, channel}.best_tperf > tol_max
            net = init(net);           
            [net_pca{i, channel}, tr_pca{i, channel}] = train(net, x, t');
            cont1 = cont1 + 1;
            if cont1 > 10 % caso não consiga atingir a meta em 5 tentativas, utilizará o ultimo treinamento realizado 
                tr_pca{i, channel}.best_tperf
                break
            end                   
         end   
    end
end
toc

disp('Treinamento completo.')
disp('...')


%% SAVE DATA
save('DADOS_TREINAMENTO\net_treinada_CIPIC.mat','net_pca', 'tr_pca');,...
%                                                 'sig', 'mu');
disp('Dados salvos!')

%% PLOT - MSE for the training set
for k = 1:no_channels
    for l = 1:no_directions
        p(l, k) = tr_pca{l, k}.best_tperf;
    end
end
   
plot(p(:,1)); hold on
plot(p(:,2)); hold off
xlim([0 1250])

legend('L', 'R');
xlabel('Indice da posição')
ylabel('RMSE')
title('Performance no mini-grupo de testes')
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
[no_subjects, no_PC, no_directions, no_channels] = size(target);


%% TREINAMENTO 
input = anthro;
trainedNet = cell(no_directions, no_channels);
nodes = 500;
in_size = size(input,1);
tic

for k = 1:1
%%% ARQUITETURA DA REDE %%%
    layers = [
        imageInputLayer([1, in_size ])          
          
        fullyConnectedLayer(nodes,'WeightsInitializer','glorot','biasInitializer', 'ones')         
        reluLayer
        
        fullyConnectedLayer(nodes,'WeightsInitializer','glorot','biasInitializer', 'ones')         
        reluLayer
   
        fullyConnectedLayer(no_PC) 
      
        regressionLayer
    ]; 
                    
    for l = 1:1    
         
        % SPLIT DATASETS for validation
        [X, PS1] = mapminmax((anthro(:, 2:end, k)), 0, 1);  
        [Y, PS2] = mapminmax(round(target(2:end, :, l, k), 2).', 0, 1);% output target varia com a direção
        
         trainRatio = 0.8;
         valRatio = 0.2;
         testRatio = 0;
        [XTraiin, XVal, ~] = divideint(X, trainRatio,valRatio,testRatio);          
        [YTraiin, YVal, ~] = divideint(Y, trainRatio,valRatio,testRatio);    
        XTrain(1,:,1,:)  = XTraiin;
        YTrain(1,1,:,:)  = YTraiin;
        XValidation(1,:,1,:)  = XVal;
        YValidation(1,1,:,:)  = YVal;       

        %%% Opções de treinamento %%%           
        options = trainingOptions('sgdm', ...
                      'MaxEpochs', 1000, ...
                      'InitialLearnRate', 0.001,... 
                      'LearnRateSchedule','piecewise',...
                      'LearnRateDropPeriod',1000,...  
                      'ValidationData', {XValidation, YValidation}, ...
                      'ValidationPatience', 10, ...
                      'Verbose', 1);                             
%                      'Plot','training-progress'); %mostra grafico
          
        %%% TREINAMENTO %%%
        [trainedNet{l,k}, trainInfo{l,k}] = trainNetwork(XTrain, YTrain, layers, options);
    end
end



%% Test network 
xinp(1, :) = mapminmax('apply', (anthro(:, 1, k)), PS1);
yout = mapminmax('reverse', predict(trainedNet{l,k}, xinp).', PS2);

figure()
trg = round(target(1, :, l, k), 2);
plot(trg); hold on;
plot(yout); hold off

legend('objetivo', 'net')
axis tight



figure()
DTF_obj = (PCWs(:, :, 1, 1) * trg.') + med_vec2(:, 1, 1); 
DTF_sim = (PCWs(:, :, 1, 1) * yout) + med_vec2(:, 1, 1); 


plot(DTF_obj); hold on;
plot(DTF_sim); hold off
axis([0 100 -50 10])
legend('objetivo', 'net')








%% SAVE DATA
path_save = 'DADOS_TREINAMENTO\deepnet_treinada';
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
if any(strcmp('riec', Datasets))
    path_save = append(path_save, '_TUB');
end
save(path_save, 'trainedNet', 'trainInfo', 'sig', 'mu');
disp('Dados Salvos!')

%% Desligar pc
% system('shutdown -s')
% disp('O sistema será desligado em 60s!')












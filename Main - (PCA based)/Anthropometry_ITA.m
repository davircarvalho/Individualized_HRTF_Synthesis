clear all; clc
%% Dados ANTHROPOMETRICOS
% Carregar dados 
path = 'B:\Documentos\#3 - TCC\EAC-TCC-Davi\HRTF-datasets\Banco ITA\'; addpath(genpath(path));
for n = 2:49
    D(:, n-1) = xlsread('Dimensions.xlsx', 1, ['I', num2str(n) ':P' num2str(n)]);
    % C: head radius, F: head depth
    x(:, n-1) = xlsread('Dimensions.xlsx', 1, ['C', num2str(n) ':D' num2str(n)]);
%     h(:, n-1) = xlsread('Dimensions.xlsx', 1, ['H', num2str(n)]);
end

% x1, x3: head width, head depth
% X = x([1,3],:)*2/10; 
X = x.*2./10;
D = D./10;
% clean inconsistencies
X(:, 14:15) = [];
D(:, 14:15) = [];

%% SALVAR DADOS
save('anthro_ITA.mat', 'X','D');
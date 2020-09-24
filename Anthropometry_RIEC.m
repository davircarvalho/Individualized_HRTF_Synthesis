clear all; clc
%% Dados ANTHROPOMETRICOS
% Carregar dados 
path = 'B:\Documentos\#3 - TCC\EAC-TCC-Davi\HRTF-datasets\Banco RIEC\'; addpath(genpath(path));
for n = 2:38
    D(:, n-1) = xlsread( 'Anthopometry_RIEC.xlsx', 1, ['B', num2str(n) ':I' num2str(n)])./10;
    X(:, n-1) = xlsread( 'Anthopometry_RIEC.xlsx', 1, ['J', num2str(n) ':L' num2str(n)])./10;
    id(n-1, :) = xlsread('Anthopometry_RIEC.xlsx', 1, ['A', num2str(n)]);
end

%% SALVAR DADOS
save('anthro_RIEC.mat', 'D', 'X', 'id');
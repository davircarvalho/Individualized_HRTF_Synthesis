clear all; clc
%% Dados ANTHROPOMETRICOS
% Carregar dados 
path = 'B:\Documentos\#3 - TCC\EAC-TCC-Davi\HRTF-datasets\Banco 3D3A lab'; addpath(genpath(path));
for n = 2:24
    D(:, n-1) = xlsread('anthro_3D3A.xls', 1, ['B', num2str(n) ':I' num2str(n)])./10;
    X(:, n-1) = xlsread('anthro_3D3A.xls', 1, ['J', num2str(n) ':L' num2str(n)])./10;
    id(n-1, :) = xlsread('anthro_3D3A.xls', 1, ['A', num2str(n)]);
end


%% SALVAR DADOS
save('anthro_3D3A.mat', 'D', 'X', 'id');
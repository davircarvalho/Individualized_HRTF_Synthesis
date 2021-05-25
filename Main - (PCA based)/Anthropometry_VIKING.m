clear all; clc
% Carregar dados 
path = '..\Datasets\VIKING\';
addpath(genpath(path));

%%
tabl = readtable('anthro_viking.xlsx');
D  = table2array(tabl(:, 4:end))./10;
X  = table2array(tabl(:, 2:3))./10;
id = table2array(tabl(:, 1));

%% SALVAR DADOS
save('anthro_VIKING.mat', 'D', 'X', 'id');
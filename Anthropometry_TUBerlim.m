clear all; clc
%% Dados ANTHROPOMETRICOS
% Carregar dados 
clc
path = 'B:\Documentos\#3 - TCC\EAC-TCC-Davi\HRTF-datasets\Banco TU Berlim\Antrhopometric measures\'; addpath(genpath(path));
Table = readtable([path 'AntrhopometricMeasures.csv']);
T = table2cell(Table);
id = cell2mat(T(:,1)); % salva id 

Tdata = str2double(strrep(T(:,2:end), ',', '.'));
%% Anthropometria
D(:, 1:8)  = Tdata(:, 14:21);
D(:, 9:16) = Tdata(:, 26:33);

theta(:, 1:2) =  Tdata(:, 24:25);
theta(:, 3:4) =  Tdata(:, 36:37);
X             =  Tdata(:, 1:13);


%% SALVAR DADOS
save('anthro_TUB.mat', 'D', 'X', 'id', 'theta');
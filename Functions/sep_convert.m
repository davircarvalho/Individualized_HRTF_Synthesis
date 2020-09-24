function [vecL, vecN] = sep_convert(vec,preci,sep,K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Função para corrigir ponto por vírgula e/ou 
%                                    criar um vetor com "k" em vez de mil.
% 
% Desenvolvido pelo professor da Engenharia Acústica 
%                                William D'Andrea Fonseca, Dr. Eng.
%
% Última atualização: 17/09/2017
%
% Entradas (vetor de números, precisão, separador, uso do "k")
%   vetor de números, exeplo: vec=600:1000:6000;
%   precisão padrão Matlab (preci): '%3.1f';
%   sepeparador: 'ponto' ou 'virgula'
%   K = 1 para o uso de 'k' no vetor de string
%
%%% Exemplos
% [L, N] = sep_convert(vec,'%3.1f','virgula',1)
% [L, ~] = sep_convert(600:1000:6000,'%3.1f','virgula',1)
% sep_convert(vec,'%4.2f','virgula')
% label = sep_convert(vec,'%4.2f')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Teste de desenvolvimento
% clear all; close all; sep = 'virgula'; preci = '%3.1f'; K = 1; vec = 600:1000:6000; vec = [{1},{2},{3}]; 

%% Input
if isempty(vec) || nargin < 1
   error('The input vector is empty, have a look.')
end
if nargin < 2; sep = 'virgula'; preci = '%3.1f'; K = 0; end
if nargin < 3; sep = 'virgula'; K = 0; end
if nargin < 4; K = 0; end
if ~isnumeric(K) || isempty(K); K = 0; end
%% Is it a cell?
if iscell(vec)
   vecMtx = cell2mat(vec);
else vecMtx = vec;
end 
%% I can only handle numbers
if ~isnumeric(vecMtx(1))
    error('I can only handle numbers, have o look to the input');
end
%% Conversion
if K == 1
    if strcmp(sep,'ponto')
        for i=1:max(size(vecMtx))
            vecN(i) = vecMtx(i)/1000;
            vecL{i} = strcat(num2str(vecN(i),preci),' k');
        end
    elseif strcmp(sep,'virgula') 
        for i=1:max(size(vecMtx))
            vecN(i) = vecMtx(i)/1000;
            vecL{i} = strcat(strrep(num2str(vecN(i),preci), '.', ','),' k');
        end
    end
else
    if strcmp(sep,'ponto')
        for i=1:max(size(vecMtx))
            vecL{i} = num2str(vecMtx(i),preci); 
        end; vecN = vecMtx;
    elseif strcmp(sep,'virgula') 
        for i=1:max(size(vecMtx))
            vecL{i} = strrep(num2str(vecMtx(i),preci), '.', ',');
        end; vecN = vecMtx;
    end
end
%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
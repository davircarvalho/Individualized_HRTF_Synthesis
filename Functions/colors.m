function cor = colors(m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Função para usar cores personalizadas
% 
% Desenvolvido pelo professor da Engenharia Acústica 
%                                William D'Andrea Fonseca, Dr. Eng.
%
% Última atualização: 04/12/2017
% 
% Compatível com Matlab R2016b
%
% m: número da cor desejada
%
% Exemplo:
% colors(1) % Para usar azul claro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cores personalizadas
    color.c0 = [0/255 0/255 0/255];        % Preto
    color.c1 = [51/255 102/255 255/255];   % Azul claro
    color.c2 = [255/255 0/255 102/255];    % Magenta claro
    color.c3 = [255/255 0/255 0/255];      % Vermelho 
    color.c4 = [0/255 30/255 255/255];     % Azul
    color.c5 = [0/255 236/255 44/255];     % Verde claro
    color.c6 = [236/255 211/255 0/255];    % Amarelo
    color.c7 = [156/255 116/255 0/255];    % Marrom
    color.c8 = [0/255 3/255 142/255];      % Azul escuro
    color.c9 = [120/255 0/255 121/255];    % Roxo
    color.c10 = [11/255 119/255 0/255];    % Verde escuro
    color.c11 = [132/255 59/255 255/255];  % Lilás
    eval(['cor = color.c' num2str(m) ';']);
    
%%% See you later.    
end
function [y_out, sig, mu] = nn_preprocess(InptMtx, varargin)
%%% Regularização de matriz entre 0 e 1 para input na ANN
% normalização da  a partir da média e da variância dos parametros entre
% diferentes exemplos, a fim deixa-los numa mesma escala 

% Exemplo 
% Processa matriz de entrada X
% [Y1, sig1, mu1] = nn_preprocess(X);

% Aplica o mesmo processamento a uma outra matriz de mesma natureza
% [Y2, ~, ~] = nn_preprocess(Y, 'apply', 'sig', sig1, 'mu', mu1);

% Recupera matriz X a sua escala original
% [X, ~, ~] = nn_preprocess(Y1, 'reverse', 'sig', sig1, 'mu', mu1);

%%% REFERENCIA %%% 
% Deep Neural Network Based HRTF Personalization Using
% Anthropometric Measurements
%% Parse Arguments
% Método de processamento
defaultMethod = 'process';
validMethods = {'process', 'apply', 'reverse'};
checkMethod = @(x) any(validatestring(x,validMethods));
% Media e desvio padrao
variancia = 'sig';
defaultVarVal = 1;
media = 'mu';
defaultMuVal = 0;

p = inputParser;
addRequired(p,'InptMtx',@isnumeric);
addOptional(p,'method',defaultMethod,checkMethod)
addParameter(p, variancia, defaultVarVal)
addParameter(p, media, defaultMuVal)
parse(p,InptMtx,varargin{:})

%% Definir media e desvio padrao
switch p.Results.method
    case 'process'
        sig = std(InptMtx, 0, 2); %desvio padrão de cada PARAMETRO entre os diferentes exemplos
        mu  = mean(InptMtx, 2);   %média de cada parâmetro
    case {'apply', 'reverse'}
        sig = p.Results.sig;
        mu = p.Results.mu;
end

%% Normalize
y_out = zeros(size(InptMtx)); % inicialização da matriz
[~, ~, no_channels] = size(InptMtx);

switch p.Results.method
% Aplicar metodo (a diferença em apply e process é apenas o valor de média e desvio padrão)
    case {'apply', 'process'}        
        for k = 1:no_channels
            y_out(:,:,k) = (1 + exp(-(InptMtx(:,:,k) - mu(:,:,k))./sig(:,:,k))).^(-1);
        end
% Voltar valores em sua escala original
    case 'reverse'
        for k = 1:no_channels
            y_out(:,:,k) = mu(:,:,k) - log(InptMtx(:,:,k).^(-1) - 1).*sig(:,:,k);
        end      
end
end
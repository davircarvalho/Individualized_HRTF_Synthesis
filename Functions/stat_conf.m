function [CI, conf] = stat_conf(input,type,alpha,mu,sigma,axistype,desv,dbref)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Função para calcular a região de confiança.
% 
% Desenvolvido pelo professor da Engenharia Acústica 
%                                William D'Andrea Fonseca, Dr. Eng.
%
% Última atualização: 04/12/2017
% 
% Compatível com Matlab R2016b
%
% input: vetor com as informações a serem estimadas.
% type: 'N' para distribuição normal e 'T' para T-Student.
% alpha: valor da confinaça, e.g., 95 (para 95%).
% mu e sigma: valores de referência para comparação com a distribuição
%             normal (geralmente mu = 0 e sigma = 1).
% axistypoe: 'nodB' para cálculo em linear
%            'dB20' para cálculo usando 20*log10(mean -+ CI*sigma/sqrt(L))
%            'dB10' para cálculo usando 10*log10(mean -+ CI*sigma/sqrt(L))
% desv: 'med' para desvio padrão da média
%       'conj' para desvio padrão do conjunto
% dbref: número de referência para usar no cálculo de dB 
%
% Exemplo:
% [Amostra.AlphaCI,conf] = stat_conf(Amostra.AlphaData,'N',80,0,1,'nodB','conj');
% [Amostra.AlphaCI,~] = stat_conf(Amostra.AlphaData,'N',90,0,1,'nodB','med');
%  Amostra.AlphaCI = stat_conf(Amostra.AlphaData,'N',90);
%  Amostra.AlphaCI = stat_conf(Amostra.AlphaData,'T',90);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Teste
% clear all; close all
% input = normrnd(0.1583975280,0.026541385465557,10,3); nargin = 1;
%% Processing

if nargin < 1
    error('I need at least the input vector.');
end
if nargin < 8 
    dbref = 1;
end
if nargin < 7 
    desv='med'; dbref = 1;
end
if nargin < 6 
    axistype = 'nodB'; desv='med';
end
if nargin < 4
    mu=0; sigma=1; axistype = 'nodB'; desv='med';
end
if nargin < 3
    alpha = 95; mu=0; sigma=1; axistype = 'nodB'; desv='med';  
end
if nargin < 2
    alpha = 95;	%  Confidence = 95%
    type = 'N'; mu=0; sigma=1; axistype = 'nodB'; desv='med';
    disp('Calcularei o desvio padrão da média, distr. Normal e CI de 95\%.')
end

%%%%% Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = size(input,1);
a = 1 - alpha/100;

if strcmp(type,'t') || strcmp(type,'T')
   CI.z = tinv([a/2,  1-a/2],L-1);	       % T-Score
   CI.calc.type = 'T-Score'; CI.calc.alpha = alpha; CI.calc.nuDF = L-1;
elseif strcmp(type,'N') || strcmp(type,'n')
   CI.z = norminv([a/2,  1-a/2],mu,sigma); % Normal
   CI.calc.type = 'Normal'; CI.calc.alpha = alpha; 
   CI.calc.muN = mu; CI.calc.sigmaN = sigma;
end

CI.confidence = alpha;
% Mean and standard deviation
CI.meanInput  = mean(input);
CI.sigmaInput = std(input);

% Confidence Intervals
if     strcmp(desv,'med')  % Desvio padrão da média
  CI.interval(:,1) = CI.z(1)*CI.sigmaInput/sqrt(L);	% Confidence Intervals 
  CI.interval(:,2) = CI.z(2)*CI.sigmaInput/sqrt(L);	% <-'
  CI.confType = 'mean';
elseif strcmp(desv,'conj') % Desvio padrão do conjunto
  CI.interval(:,1) = CI.z(1)*CI.sigmaInput;	% Confidence Intervals 
  CI.interval(:,2) = CI.z(2)*CI.sigmaInput;	% <-'
  CI.confType = 'sample';  
end

% Mean +/- Standard Error
if strcmp(axistype,'lin') || strcmp(axistype,'nodb') || strcmp(axistype,'nodB')
for i=1:length(CI.meanInput)
  CI.muCI(i,1) = CI.meanInput(1,i) - abs(CI.interval(i,1));
  CI.muCI(i,2) = CI.meanInput(1,i) + abs(CI.interval(i,2));
end
%%%% decibel stuff
elseif strcmp(axistype,'dB20') || strcmp(axistype,'db20') || strcmp(axistype,'DB20') 
    CI.meanInputdB = 20.*log10(mean(input)./dbref);
    % Use 20*log10(mean -+ CI*sigma/sqrt(L))
    for i=1:length(CI.meanInput)
      CI.muCI(i,2) = 20.*log10((CI.meanInput(1,i) + abs(CI.interval(i,2)))./dbref);
      CI.muCI(i,1) = CI.meanInputdB(i) - (CI.muCI(i,2) - CI.meanInputdB(i)'); 
    end
    CI.intervaldB(:,2) = CI.muCI(:,2) - CI.meanInputdB';
    CI.intervaldB(:,1) = -CI.intervaldB(:,2);    
elseif strcmp(axistype,'dB10') || strcmp(axistype,'db10') || strcmp(axistype,'DB10') 
    CI.meanInputdB = 10.*log10(mean(input)./dbref);
    % Use 10*log10(mean -+ CI*sigma/sqrt(L))    
    for i=1:length(CI.meanInput)
      CI.muCI(i,2) = 10.*log10((CI.meanInput(1,i) + abs(CI.interval(i,2)))./dbref);
      CI.muCI(i,1) = CI.meanInputdB(i) - (CI.muCI(i,2) - CI.meanInputdB(i)'); 
    end
    CI.intervaldB(:,2) = CI.muCI(:,2) - CI.meanInputdB';
    CI.intervaldB(:,1) = -CI.intervaldB(:,2);    
end  
    
conf = abs(CI.z(2));

end
function arruma_fig(preciX,preciY,separador,figura,KX,KY,Getlabel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Função para arrumar os eixos da figura.
% 
% Desenvolvido pelo professor da Engenharia Acústica 
%                                William D'Andrea Fonseca, Dr. Eng.
%
% Última atualização: 17/09/2017
%
% Entradas (precisãoX, precisãoY, separador, figura)
%   preciX, preciY = '%3.1f' (padrão Matlab).
%   sepeparador = 'ponto' ou 'virgula' (opcional).
%   figura = objeto da figura a sert ajustada (opcional).
%   KX e KY = 0 ou 1 para escrever o numero de forma compacta com um 'k' no
%   final (opcional).
%
%%% Exemplos:
% arruma_fig('% 4.0f','% 2.2f','virgula',fig.Impedancia,1,0,[1,0])
% arruma_fig('% 4.0f','% 2.2f','virgula',fig.Impedancia,1,0)
% arruma_fig('% 4.0f','% 2.2f','virgula',fig.Impedancia)
% arruma_fig('% 4.0f','% 2.2f','virgula',gcf)
% arruma_fig('% 4.0f','% 2.2f','virgula')
% arruma_fig('% 4.0f','% 2.2f')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Development test
% preciX = '% 4.0f'; preciY = '% 2.2f'; separador = 'virgula'; figura = gcf; KX = 0; KY = 0;

%% Corrige entradas
if nargin<2; error('Declare as precisões dos eixos.'); end;
if nargin<3; separador = 'virgula'; figura = gcf; ax = figura.CurrentAxes; KX = 0; KY = 0; end; 
if ~exist('separador','var') || isempty(separador); separador = 'virgula'; end;
if nargin<4; figura = gcf; ax = figura.CurrentAxes; KX = 0; KY = 0; Getlabel= [0, 0]; end; 
if nargin<7; figura = gcf; ax = figura.CurrentAxes; KX = 0; KY = 0; Getlabel= [0, 0]; end; 
if ~exist('figura','var') || isempty(figura); figura = gcf; ax = figura.CurrentAxes; end;
if ~exist('ax','var') || isempty(ax); ax = figura.CurrentAxes; end;

%% Arrumando

%%% Eixo X
xtick.GetNum = ax.XTick; xtick.GetLabel = ax.XTickLabel;
if Getlabel(1) == 1
   xtick.NumFromLabel = str2double(xtick.GetLabel);
  [xtick.NewLabel,~]=sep_convert(xtick.NumFromLabel,preciX,separador,KX);
else
   [xtick.NewLabel,~]=sep_convert(xtick.GetNum,preciX,separador,KX);
end

set(ax,'XTick',xtick.GetNum); set(ax,'XTickLabel',xtick.NewLabel);

%%% Eixo Y
ytick.GetNum = ax.YTick; ytick.GetLabel = ax.YTickLabel; 
if Getlabel(2) == 1
   ytick.NumFromLabel = str2double(ytick.GetLabel);
   [ytick.NewLabel,~]=sep_convert(ytick.NumFromLabel,preciY,separador,KY);
else
  [ytick.NewLabel,~]=sep_convert(ytick.GetNum,preciY,separador,KY); 
end

 set(ax,'YTick',ytick.GetNum); set(ax,'YTickLabel',ytick.NewLabel);
%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
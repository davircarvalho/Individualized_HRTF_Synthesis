function [Obj_out, error] = sofaFit2Grid(Obj_in, out_pos, varargin)
% Converte posi��es de HRIRs SOFA �s posi��es especificadas em 'out_pos'
% Davi R. Carvalho @UFSM - Engenharia Acustica - julho/2020

%  ~Input Parameters:
%    Obj_in:     Objeto de HRTFs SOFA com coordenadas esf�ricas 
%                azi:    0� -> 360�     
%                elev: -90� -> 90� 
%                (com o �ltimo update, talvez funcione com outros sistemas de coordenadas)
%    out_pos:    Nx3 matrix de posi��es desejadas, em que N corresponde ao 
%                n�mero total de posi��es, e as colunas correspondem a azimute,
%                eleva��o e raio respectivamente.

%  ~Output Parameters:
%     Obj_out:   Objeto de HRTFs SOFA com as caracter�stica 
%                de medi��o do dataset CIPIC.
%     error:     Erro rms (�) para cara posicao (valido apenas em 'adapt') 

%  ~Optional Parameters:     
%    'adapt':    Seleciona as posicoes mais proximas do grid objetivo e
%                for�a a assumirem suas coordenadas.
%    'hybrid':   Faz a adapta��o como em 'adapt', mas posi��es do grid
%                original escolhidas para mais de uma posi��o objetivo 
%                s�o determinadas por interpola��o (Metodo Padrao).
%    'vbap':     Interpola��o vbap
%    'bilinear': Interpola��o bilinear 
%    'spherical_harmonics': Interpola��o por meio de harmonicos esf�ricos
%    'Fs':       Transforma��o da taxa de amostragem no objeto de sa�da 
%                (Padrao: original do objeto).

% Exemplo: Obj_out = sofaFit2Grid(Obj_in, out_pos, 'spherical_harmonics', 'Fs', 48000)
%
% Matlab R2020a
%% Parse Arguments
% M�todo de processamento
defaultMethod = 'adapt';
validMethods = {'adapt', 'hybrid', 'vbap', 'bilinear', 'spherical_harmonics', 'sh', ...
    'hybrid2'};
checkMethod = @(x) any(validatestring(x,validMethods));

% Op��es de taxa de amostragem
paramName = 'Fs';
defaultVal = Obj_in.Data.SamplingRate;

%Verificar entradas
p = inputParser;
addRequired(p,'Obj_in',@isstruct);
addOptional(p,'method',defaultMethod,checkMethod)
addParameter(p,paramName,defaultVal)
parse(p,Obj_in,varargin{:})


%% Sample rate match
if Obj_in.Data.SamplingRate ~= p.Results.Fs
    Obj_in = sofaResample(Obj_in, p.Results.Fs);
end

%% Initialize error
error = zeros(length(out_pos),1);
meta.pos = Obj_in.SourcePosition;
idx_adapt = zeros(length(out_pos), 1);
for zz = 1:length(out_pos) 
    % Calculo do erro entre posi��o objetivo e posi��es disponiveis
    tsqr = sqrt((meta.pos(:,1)-out_pos(zz,1)).^2 + (meta.pos(:,2)-out_pos(zz,2)).^2);
    [error(zz), idx_adapt(zz,1)] = min(tsqr); 
    idx_adapt(zz,2) = zz; %salvar indice da posi��o objetivo          
end
meta.fittedIR = Obj_in.Data.IR(idx_adapt(:,1), :, :);

idx_hybrid = find(error > 0.5); % find indexes where the error is "high"

% thr = 55;
% idx_hybrid = find(out_pos(:, 2)<thr & out_pos(:, 2)>-thr);
des_hybrid = out_pos(idx_hybrid,:);
IR_temp = interpolateHRTF(Obj_in.Data.IR, meta.pos(:,[1,2]), des_hybrid(:,[1,2]), ...
                          'Algorithm','bilinear');
meta.fittedIR(idx_hybrid,:,:) = IR_temp;

%% OUTPUT data (assembly and metadata) 
Obj_out = SOFAgetConventions('SimpleFreeFieldHRIR');
Obj_out.Data.IR = meta.fittedIR;
Obj_out.SourcePosition = out_pos;
Obj_out.Data.SamplingRate = p.Results.Fs;

% warning('off','SOFA:upgrade');
% Obj_out = SOFAupgradeConventions(Obj_out);
Obj_out = SOFAupdateDimensions(Obj_out);  
end

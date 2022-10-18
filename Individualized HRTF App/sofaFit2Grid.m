function [Obj_out, error] = sofaFit2Grid(Obj_in, out_pos, varargin)
% Converte posições de HRIRs SOFA às posições especificadas em 'out_pos'
% Davi R. Carvalho @UFSM - Engenharia Acustica - julho/2020

%  ~Input Parameters:
%    Obj_in:     Objeto de HRTFs SOFA com coordenadas esféricas 
%                azi:    0° -> 360°     
%                elev: -90° -> 90° 
%                (com o último update, talvez funcione com outros sistemas de coordenadas)
%    out_pos:    Nx3 matrix de posições desejadas, em que N corresponde ao 
%                número total de posições, e as colunas correspondem a azimute,
%                elevação e raio respectivamente.

%  ~Output Parameters:
%     Obj_out:   Objeto de HRTFs SOFA com as característica 
%                de medição do dataset CIPIC.
%     error:     Erro rms (°) para cara posicao (valido apenas em 'adapt') 

%  ~Optional Parameters:     
%    'adapt':    Seleciona as posicoes mais proximas do grid objetivo e
%                força a assumirem suas coordenadas.
%    'hybrid':   Faz a adaptação como em 'adapt', mas posições do grid
%                original escolhidas para mais de uma posição objetivo 
%                são determinadas por interpolação (Metodo Padrao).
%    'vbap':     Interpolação vbap
%    'bilinear': Interpolação bilinear 
%    'spherical_harmonics': Interpolação por meio de harmonicos esféricos
%    'Fs':       Transformação da taxa de amostragem no objeto de saída 
%                (Padrao: original do objeto).

% Exemplo: Obj_out = sofaFit2Grid(Obj_in, out_pos, 'spherical_harmonics', 'Fs', 48000)
%
% Matlab R2020a
%% Parse Arguments
% Método de processamento
defaultMethod = 'adapt';
validMethods = {'adapt', 'hybrid', 'vbap', 'bilinear', 'spherical_harmonics', 'sh', ...
    'hybrid2'};
checkMethod = @(x) any(validatestring(x,validMethods));

% Opções de taxa de amostragem
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
    % Calculo do erro entre posição objetivo e posições disponiveis
    tsqr = sqrt((meta.pos(:,1)-out_pos(zz,1)).^2 + (meta.pos(:,2)-out_pos(zz,2)).^2);
    [error(zz), idx_adapt(zz,1)] = min(tsqr); 
    idx_adapt(zz,2) = zz; %salvar indice da posição objetivo          
end
meta.fittedIR = Obj_in.Data.IR(idx_adapt(:,1), :, :);

thr = 55;
idx_hybrid = find(out_pos(:, 2)<thr & out_pos(:, 2)>-thr);
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

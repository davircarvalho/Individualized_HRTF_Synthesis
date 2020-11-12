function Obj_out = sofaFit2Grid(Obj_in, out_pos, varargin)
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

% Exemplo: Obj_out = sofaFit2Grid(Obj_in, out_pos, 'bilinear', 'Fs', 48000)
%
% Matlab R2020a
%% Parse Arguments
% Método de processamento
defaultMethod = 'adapt';
validMethods = {'adapt','hybrid', 'vbap', 'bilinear', 'spherical_harmonics'};
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


%% "Interpolation" by nearest position ('ADAPT')
switch p.Results.method
    case {validMethods{1}, validMethods{2}}              
        meta.pos = Obj_in.SourcePosition;
        idx_adapt = zeros(length(out_pos), 1);
        meta.fittedPOS = zeros(size(out_pos));
        for zz = 1:length(out_pos) 
            % Calculo do erro entre posição objetivo e posições disponiveis
            tsqr = sqrt((meta.pos(:,1)-out_pos(zz,1)).^2 + (meta.pos(:,2)-out_pos(zz,2)).^2);
            [~, idx_adapt(zz,1)] = min(tsqr); 
            idx_adapt(zz,2) = zz; %salvar indice da posição objetivo
            
            % Posicoes selecionadas no grid original (util para visualização)
%             meta.fittedPOS(zz,:) = Obj_in.SourcePosition(idx_adapt(zz, 1),:);            
        end
        % (caso for visualizar as posições selecionadas, comentar a linha abaixo)
        meta.fittedPOS = out_pos; % <-------
        meta.fittedIR  = Obj_in.Data.IR(idx_adapt(:,1), :, :);
        
        %% Modelo Hibrido ('HYBRID') 
        if any(strcmp(validMethods{2}, p.Results.method)) 
            % selecionar apenas valores repetidos 
            idx_hybrid = idx_adapt;
            [~,ind_uniq] = unique(idx_adapt(:,1));
            idx_hybrid = removerows(idx_hybrid, 'ind', ind_uniq);
             if isempty(idx_hybrid) %sem indice repetido, sem interpolacao
        %         warning('Nenhum índice repetido identificado, método apenas adaptativo.')
             else         
                meta.fittedIR(idx_hybrid(:,2),:,:) = NaN; % to be filled later           
                % interpolar valores limpos
                des_hybrid = out_pos(idx_hybrid(:,2),:);
                IR_temp = interpolateHRTF(Obj_in.Data.IR, meta.pos(:,[1,2]), des_hybrid(:,[1,2]), ...
                                          'Algorithm','bilinear');   
                meta.fittedIR(idx_hybrid(:,2),:,:) = IR_temp;
             end
        end

%% Interpolar por 'VBAP' ou 'BILINEAR'
    case {validMethods{3}, validMethods{4}}
        meta.fittedIR = zeros(length(out_pos), 2, size(Obj_in.Data,3));
        meta.pos = Obj_in.SourcePosition;
        meta.fittedIR = miinterpolateHRTF(Obj_in.Data.IR, meta.pos(:,[1,2]), out_pos(:,[1,2]),...
                                       'Algorithm', p.Results.method);   
        meta.fittedPOS = out_pos;

%% Interpolar por harmonicos esféricos 'spherical_harmonics'
    case {validMethods{5}}
        IR_interp = sofaSHinterpolate(Obj_in, out_pos(:, [1,2]));
        meta.fittedIR = IR_interp.Data.IR;
        meta.fittedPOS = out_pos;
end 



%% OUTPUT data (assembly and metadata) 
Obj_out = SOFAgetConventions('SimpleFreeFieldHRIR');
Obj_out.Data.IR = meta.fittedIR;
Obj_out.SourcePosition = meta.fittedPOS;
Obj_out.Data.SamplingRate = p.Results.Fs;

% warning('off','SOFA:upgrade');
Obj_out = SOFAupgradeConventions(Obj_out);
Obj_out = SOFAupdateDimensions(Obj_out);


%% Plots
% ri = Obj_in.SourcePosition(1,3);
% %%% plot input %%%
% Obj_in.SourceView = Obj_in.SourcePosition;
% Obj_in.SourceView_Type = 'spherical';
% Obj_in.API.Dimensions.SourceView  = 'MC';
% SOFAplotGeometry(Obj_in)
% 
% view([35 20])
% % xlabel('X [m]')
% % ylabel('Y [m]')
% % zlabel('Z [m]')
% % xticks([-ri, 0, ri])
% % yticks([-ri, 0, ri])
% % zticks([-ri, 0, ri])
% % xticklabels([-ri, 0, ri])
% % yticklabels([-ri, 0, ri])
% % zticklabels([-ri, 0, ri])
% set(gca,'XColor', 'none','YColor','none', 'ZColor','none')
% name = 'CIPIC';
% title(name)
% legend off
% axis tight
% export_fig([pwd, '\Images\English\' name ], '-pdf', '-transparent');

% % %%% plot output %%%ri = Obj_in.SourcePosition(1,3);
% ro = Obj_in.SourcePosition(1,3);
% Obj_out.SourceView = Obj_out.SourcePosition;
% Obj_out.SourceView_Type = 'spherical';
% Obj_out.API.Dimensions.SourceView  = 'MC';
% SOFAplotGeometry(Obj_out)
% view([35 20])
% xlabel('X [m]')
% ylabel('Y [m]')
% zlabel('Z [m]')
% xticks([-ro, 0, ro])
% yticks([-ro, 0, ro])
% zticks([-ro, 0, ro])
% xticklabels([-ro, 0, ro])
% yticklabels([-ro, 0, ro])
% zticklabels([-ro, 0, ro])
% title('')
% legend off
% % axis tight
% % export_fig([pwd, '\Images\3dITAout' ], '-pdf', '-transparent');


         
% % plot error (Mapa 2d sobreposto) %--------------------------------------
% in_pos = Obj_in.SourcePosition(idx_adapt(:,1), :);
% figure()
% scatter(out_pos(:,1), out_pos(:,2), 27,'k', 'filled'); hold on 
% scatter(in_pos(:,1), in_pos(:,2), 25, 'r', 'filled', 'square'); hold off
% xlabel('Azimute [grau]')
% ylabel('Elevação [grau]')
% legend('Objetivo', 'Original', 'location', 'southeast')
% axis tight
% set(gca,'FontSize',12)
        

end

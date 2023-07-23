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
validMethods = {'adapt', 'hybrid', 'vbap', 'bilinear', 'spherical_harmonics', 'sh'};
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

%% "Interpolation" by nearest position ('ADAPT')
switch p.Results.method
    case {validMethods{1}, validMethods{2}}              
        meta.pos = Obj_in.SourcePosition;
        idx_adapt = zeros(length(out_pos), 1);
        meta.fittedPOS = zeros(size(out_pos));
        for zz = 1:length(out_pos) 
            % Calculo do erro entre posi��o objetivo e posi��es disponiveis
            tsqr = sqrt((meta.pos(:,1)-out_pos(zz,1)).^2 + (meta.pos(:,2)-out_pos(zz,2)).^2);
            [error(zz), idx_adapt(zz,1)] = min(tsqr); 
%             idx_adapt(zz,1) = dsearchn(meta.pos(:,[1,2]), out_pos(zz,[1,2]));
            idx_adapt(zz,2) = zz; %salvar indice da posi��o objetivo
            
            % Posicoes selecionadas no grid original (util para visualiza��o)
%             meta.fittedPOS(zz,:) = Obj_in.SourcePosition(idx_adapt(zz, 1),:);            
        end
        meta.fittedIR  = Obj_in.Data.IR(idx_adapt(:,1), :, :);
        
        %% Modelo Hibrido ('HYBRID') 
        if any(strcmp(validMethods{2}, p.Results.method)) 
            % selecionar apenas valores repetidos 
%             idx_hybrid = idx_adapt;
%             [~,ind_uniq] = unique(idx_adapt(:,1));
             idx_hybrid = find(error > 1);
             size(idx_hybrid)
             if isempty(idx_hybrid) %sem indice repetido, sem interpolacao
        %         warning('Nenhum �ndice repetido identificado, m�todo apenas adaptativo.')
             else         
%                 meta.fittedIR(idx_hybrid,:,:) = NaN; % to be filled later           
                % interpolar valores limpos
                des_hybrid = out_pos(idx_hybrid,:);
                try
                    IR_temp = interpolateHRTF(Obj_in.Data.IR, meta.pos(:,[1,2]), des_hybrid(:,[1,2]), ...
                                              'Algorithm','bilinear');   
                    meta.fittedIR(idx_hybrid,:,:) = IR_temp;
                catch
                end
             end
        end

%% Interpolar por 'VBAP' ou 'BILINEAR'
    case {validMethods{3}, validMethods{4}}
        meta.fittedIR = zeros(length(out_pos), 2, size(Obj_in.Data,3));
        meta.pos = Obj_in.SourcePosition;
        meta.fittedIR = interpolateHRTF(Obj_in.Data.IR, meta.pos(:,[1,2]), out_pos(:,[1,2]),...
                                          'Algorithm', p.Results.method);   

%% Interpolar por harmonicos esf�ricos 'spherical_harmonics'
    case {validMethods{5}, validMethods{6}}
        IR_interp = sofaSHinterpolate(Obj_in, out_pos(:, [1,2]));
        meta.fittedIR = IR_interp.Data.IR;
end 


%% OUTPUT data (assembly and metadata) 
Obj_out = SOFAgetConventions('SimpleFreeFieldHRIR');
Obj_out.Data.IR = meta.fittedIR;
Obj_out.SourcePosition = out_pos;
Obj_out.Data.SamplingRate = p.Results.Fs;

% warning('off','SOFA:upgrade');
% Obj_out = SOFAupgradeConventions(Obj_out);
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
% ylabel('Eleva��o [grau]')
% legend('Objetivo', 'Original', 'location', 'southeast')
% axis tight
% set(gca,'FontSize',12)
        

end

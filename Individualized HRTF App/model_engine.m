% %for test purposes only
% x1 = 17;
% x3 = 20;
% d1 = 1;
% d3 = 3;
% d5 = 5;
% d6 = 4;
% d7 = 1;
% d8 = 1;
% FileName = 'EAC_test';
% hesuvi = 1;
% sofa = 0;
% Fs_out = 96000;
% FR = {'30°, 0°'};
% SR = {'90°, -30°'};
% RR = {'150°, 0°'};
% FF = {'0°, 15°'};
% FL = {'30°, 0°'};
% SL = {'90°, 0°'};
% RL = {'160°, -15°'};

% DAVI ROCHA CARVALHO; ENG. ACUSTICA - UFSM; AGOSTO/2020
% Simulação de HRTF individualizada, a partir da RNA treinada e uso de dados
% antropométricos
function model_engine(x1, x3,... % anthro cabeca
                      d1L, d2L, d3L, d5L, d7L, d8L, ...   % antropometria L
                      d1R, d2R, d3R, d5R, d7R, d8R, ...   % antropometria R
                      FileName, hesuvi_flag, sofa_flag, Fs_out, ... % useful inputs
                      FR, SR, RR, FF, FL, SL, RL)         %HeSuVi angles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load networks
tic
load('net_treinada_CIPIC_ARI_ITA_3D3A.mat');
load('target_pca_CIPIC_ARI_ITA_3D3A.mat');
Fs_sim = 44100; % sample rate das IR geradas na nn  

%% Defnindo INPUT
% separando os 8 parametros necessários para a ann
dL = [d1L, d2L, d3L, d5L, d7L, d8L];
dR = [d1R, d2R, d3R, d5R, d7R, d8R];

x = [x1, x3];

% unir as duas matrizes
InptMtx(:,:,1) = abs([x, dL]');  % L 
InptMtx(:,:,2) = abs([x, dR]');  % R

[~, ~, no_directions, no_channels] = size(PCWs);

%% Simulação do modelo
wait = waitbar(0,'Neural net inference ..');
cont = 0;
for n = 1:no_channels
    for i = 1:no_directions
        result = net_pca{i, n}(InptMtx(:,:, n));
        DTF_sim(:, i, n) = (PCWs(:, :, i, n) * result) + med_vec2(:, i, n); 
        
        cont = cont+1;
        waitbar((cont)/(no_directions*no_channels), wait);        
    end
end

waitbar(0.3, wait, 'Processing results ...');

%% RECONSTRUÇÃO DE FASE: (fase mínima + ITD)
% CÁlCULO DO ITD
new_itd = itd_synthesis(x(1), x(2), out_pos, Fs_sim, 'adapt');
for k = 1:no_directions
    itd = new_itd(k);
    offset = 20.4;
    [IR_minL, IR_minR] = phase_job(DTF_sim(:, k, 1), DTF_sim(:, k, 2), ...
                                   itd, out_pos(k,:), offset); 
    %save no formato CIPIC
    hrir_final(:, k, 1) = IR_minL;
    hrir_final(:, k, 2) = IR_minR;
end

%% EXPORTS
Obj = SOFAgetConventions('SimpleFreeFieldHRIR');
Obj.Data.IR = shiftdim(hrir_final, 1);
Obj.Data.SamplingRate = Fs_sim; % valor atual 
Obj.SourcePosition = out_pos;
Obj = SOFAupdateDimensions(Obj);

%%% save HeSuVi %%%%%%%
if hesuvi_flag 
    waitbar(0.8,wait, 'Assembling HeSuVi output .....');
    hesuvi_angles = {FR; SR; RR; FF; FL; SL; RL};
    HeSuViExport(Obj, FileName, Fs_out, hesuvi_angles)
end

%%% save SOFA file %%%%
if sofa_flag 
    waitbar(0.6,wait,'Assembling SOFA output ....');
    if Fs_out ~= Fs_sim 
        Obj = sofaResample(Obj, Fs_out);
    end
    file_path = ([FileName, '_', num2str(Fs_out/1000), 'kHz.sofa']);
    
    compression = 0;
    SOFAsave(file_path, Obj, compression);
end


waitbar(1, wait, ['All Done!          Time: ' num2str(toc) 's']);

pause(2)
close(wait)
end











%% INTERNAL FUNCTIONS
function HeSuViExport(Obj_in, FileName, Fs_out, angles)
    %% Find index of speaker angles
    for k = 1:length(angles)
        A = regexp(angles{k}, '[+-]?\d+', 'match');
         temp = cellfun(@str2num, A);
         azim(k) = temp(1);
         elev(k) = temp(2);
        if k >= 5
            azim(k) = 360 - azim(k);
        end 
    end
    
    
%%% Assemble and Processing %%% -------------------------------------------
    % Get index of measurements with the same directions   
    Obj_in.SourcePosition = round(Obj_in.SourcePosition);
    n = 0;
    for k = 1:length(azim)
        tsqr = sqrt((Obj_in.SourcePosition(:,1) - azim(k)).^2 + (Obj_in.SourcePosition(:,2) - elev(k)).^2);
        [~,pos_idx] = min(tsqr); 
%         pos_idx = find(Obj_in.SourcePosition(:,1)==azim(k) & Obj_in.SourcePosition(:,2)==elev(k)); %#ok<SAGROW>
        
        clear Obj_out IR
        Obj_out = SOFAgetConventions('SimpleFreeFieldHRIR'); 
        Obj_out.Data.SamplingRate = Obj_in.Data.SamplingRate;
        Obj_out.Data.IR = Obj_in.Data.IR(pos_idx,:,:);
        Obj_out.SourcePosition = Obj_in.SourcePosition(pos_idx,:);
        Obj_out = SOFAupdateDimensions(Obj_out);
        
        % processing IR 
        if Fs_out ~= Obj_in.Data.SamplingRate
            Obj_out = sofaResample(Obj_out, Fs_out, 4096); %tamanho do vetor de saida pra permitir low freq
        end
        Obj_out = sofaAddLowFrequency(Obj_out);

        % prepare for hesuvi export
        IR = squeeze(Obj_out.Data.IR).'; 
        for g = 1:size(IR, 2) % L/R channels
            n = n+1;
            y(:,n) = IR(:,g);
        end
    end
%%% Export to .wav %%%-----------------------------------------------------
    y = y./max(abs(y(:)));
    y = y(:, [1,2,3,4,5,6,7,10,9,12,11,14,13,8]);
    file_path = ([FileName, '_', num2str(Fs_out/1000), 'kHz.wav']);
    audiowrite(file_path, y, Obj_out.Data.SamplingRate)
end


%% INTERNAL FUNCTIONS
function Obj = sofaResample(Obj, Fs, Nintp)
    N = ceil( (Fs/Obj.Data.SamplingRate) * size(Obj.Data.IR, 3) ); % length after resample
    if nargin < 3  ||  Nintp < N
        Nintp = 2^nextpow2(N); % output length
    end
    zpad = zeros((Nintp - N), 1);
    % options
    [p,q] = rat(Fs / Obj.Data.SamplingRate);
    normFc = .98 / max(p,q);
    order = 256 * max(p,q);
    beta = 12;
    %%% Cria um filtro via Least-square linear-phase FIR filter design
    lpFilt = firls(order, [0 normFc normFc 1],[1 1 0 0]);
    lpFilt = lpFilt .* kaiser(order+1,beta)';
    lpFilt = lpFilt / sum(lpFilt);
    % multiply by p
    lpFilt = p * lpFilt;
    % Actual Resample
    for k = 1:size(Obj.Data.IR, 1)
        for l = 1:size(Obj.Data.IR, 2)
            IRpre(k, l, :) = resample(Obj.Data.IR(k, l, :),p,q,lpFilt);
            IR(k, l, :) = [squeeze(IRpre(k, l, :)); zpad];
        end 
    end
    Obj.Data.IR = IR;
    Obj.Data.SamplingRate = Fs;
end
%%


%   %Trocar mensagem do botão
%             app.IniciarButton.Text = 'Processando..';
%                 
%             % check requested inputs
%             sofa_flag = app.ExportSOFACheckBox.Value;
%             hesuvi_flag = app.ExportHeSuVi71wavCheckBox.Value;
%             
%             %HeSuVi angles
%             % FR, SR, RR, RL, SL, FL, FF -> current order
%             % FR, SR, RR, FF, FL, SL, RL  -> objective order
%            
%             % Run Engine
%             model_engine(app.x1.Value, app.x3.Value, app.d1.Value, ...
%                 app.d3.Value, app.d5.Value, app.d6.Value, ... 
%                 app.d7.Value, app.d8.Value, app.ID.Value, ...
%                 hesuvi_flag, sofa_flag, ...
%                 app.SampleRateHzopcional.Value,...
%                 app.AzimElevDropDown.Value,...    
%                 app.AzimElevDropDown_2.Value,... 
%                 app.AzimElevDropDown_3.Value,... 
%                 app.AzimElevDropDown_7.Value,...
%                 app.AzimElevDropDown_6.Value,... 
%                 app.AzimElevDropDown_5.Value,...
%                 app.AzimElevDropDown_4.Value)
%                
%             app.IniciarButton.Text = 'Iniciar';

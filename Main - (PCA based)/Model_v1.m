clear all; close all; clc
% DAVI ROCHA CARVALHO; ENG. ACUSTICA - UFSM; JUNHO/2019
% Simulação de HRTF individualizada, a partir da RNA treinada e uso de dados
% antropom�tricos
addpath([pwd, '\..\Functions']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Options
% Defina quais datasets usar: {'cipic', 'ari', 'ita', '3d3a', 'riec'}
Datasets = {'cipic', 'ari', 'ita', '3d3a'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load networks
net_file = 'DADOS_TREINAMENTO/net_treinada';
target_file = 'DADOS_TREINAMENTO/target_pca';
if any(strcmp('cipic', Datasets))
    net_file = append(net_file, '_CIPIC');
    target_file = append(target_file, '_CIPIC');
end
if any(strcmp({'ari'}, Datasets))
    net_file = append(net_file, '_ARI');
    target_file = append(target_file, '_ARI');
end
if any(strcmp({'ita'}, Datasets))
    net_file = append(net_file, '_ITA');
    target_file = append(target_file, '_ITA');
end
if any(strcmp({'3d3a'}, Datasets))
    net_file = append(net_file, '_3D3A');
    target_file = append(target_file, '_3D3A');
end
if any(strcmp({'riec'}, Datasets))
    net_file = append(net_file, '_RIEC');
    target_file = append(target_file, '_RIEC');
end
if any(strcmp({'tub_meas'}, Datasets))
    net_file = append(net_file, '_TUBMEAS');
    target_file = append(target_file, '_TUBMEAS');
end
if any(strcmp({'tub_sim'}, Datasets))
    net_file = append(net_file, '_TUBSIM');
    target_file = append(target_file, '_TUBSIM');
end

load(net_file)
load(target_file)
fs = 44100;

disp('Network carregada!')

%% Defnindo INPUT --------------------------------------------------------------
% ID
subject = 140;
subj = num2str(subject);

% separando os 8 parametros necess�rios para a ann
% d1:d8
dL = [1.98, 0.845, 1.75, 2.02, 7.42, 4.185, 0.68, 1.3];
d1 = dL([1,2,3,5,7,8]);
d2 = d1;

% x1, x3, x12
% x1 = head width
% x3 = head depth
% x12 = shoulder width
% x = anthro.X(subject, [1, 3, 12]);
x = [16.3, 19.5, 46];

% unir as duas matrizes
xinpt(:,:,1) = abs([x([1,2]), d1]');  % L
xinpt(:,:,2) = abs([x([1,2]), d2]');  % R

[no_samples, no_PC, no_directions, no_channels] = size(PCWs);

%% Simulação do modelo
disp('...')
disp('Simulação iniciada')
% computar valores de saida para dada entrada na ann
tic
wait = waitbar(0,'Processando Novas Componentes Principais');
cont = 0;
for n = 1:no_channels
    for i = 1:no_directions
        result = net_pca{i, n}(xinpt(:,:, n));
        DTF_sim(:, i, n) = (PCWs(:, :, i, n) * result) + med_vec2(:, i, n);        
        cont = cont+1;
        waitbar((cont)/(no_directions*no_channels), wait);        
    end
end
toc
close(wait)
disp('Simulação concluida!')
disp('...')

%% RECONSTRUÇÃO DE FASE: (fase mínima + ITD)
% CÁlCULO DO ITD
new_itd = itd_synthesis(x(1), x(2), out_pos, fs, 'adapt');

% Aplicando fase minima e de excesso 
for k = 1:no_directions
    itd = new_itd(k);
    offset = 20;  
    [IR_minL, IR_minR] = phase_job(DTF_sim(:, k, 1), DTF_sim(:, k, 2), ...
                                   itd, out_pos(k,:), offset); 
    %save no formato CIPIC
    hrir_final(:, k, 1) = IR_minL;
    hrir_final(:, k, 2) = IR_minR;
end


%% Export para SOFA
Obj = SOFAgetConventions('SimpleFreeFieldHRIR');
Obj.Data.IR = shiftdim(hrir_final, 1);
Obj.Data.SamplingRate = fs;
Obj.SourcePosition = out_pos;
Obj = SOFAupdateDimensions(Obj);

%%%% save SOFA file %%%%
    file_path = fullfile(['individuo_' subj '.sofa']);
    disp(['Saving:  ' file_path])
    compression = 0;
    SOFAsave(file_path, Obj, compression);

%%%% save HeSuVi %%%%
HeSuViExport(Obj)

disp('all done!')
toc


%% Plot geometry mesh
% Obj.SourceView_Type = 'spherical';
% Obj.API.Dimensions.SourceView  = 'MC';
% Obj.SourceView = Obj.SourcePosition;
% SOFAplotGeometry(Obj)


%% PLOT - HRIR [Simulada vs Medida]
% figure()
%
% % DEFINA A DIRE��O DA RI
% azim = 90;
% elev = 0;
%
% % Get index of measurements with the same directions
% pos=find(round(Obj.SourcePosition(:,1))==azim & round(Obj.SourcePosition(:,2))==elev);
%
%
% %%% ESQUERDA --------------------------------------------------------------
% subplot(2,1,1)
% plot(real(squeeze(Obj.Data.IR(pos, 1, :))), 'linewidth', 1.5);
%
% xlabel('Amostras')
% ylabel('Esquerda')
% azim = num2str(out_pos(pos,1));
% elee = num2str(out_pos(pos,2));
% title(['HRIR Invididuo ', num2str(subject), ', azimute ', azim, ...
%                                    '�, eleva��o ', elee, '�']);
% set(gca,'FontSize',12)
% axis tight
%
% %%% DIREITA ---------------------------------------------------------------
% subplot(2,1,2)
% plot( real(squeeze(Obj.Data.IR(pos, 2, :))), 'linewidth', 1.5); hold on
%
% xlabel('Amostras')
% ylabel('Direita')
%
% % ylim([-1 1])
% set(gca,'FontSize',12)
% axis tight


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT - HRTF_log para determinada dire��o �
% N = no_samples;
% freq = linspace(0, fs-fs/N, N);
% f_1000hz = dsearchn(freq', 1000); %posicao da sample em 1000Hz
%
% %%%% DEFINA A DIRE��O DA HRTF %%%
% azim = 90;
% elev = 0;
%
% % Get index of measurements with the same directions
% pos=find(Obj.SourcePosition(:,1)==azim & Obj.SourcePosition(:,2)==elev);
%
% %%% Simulado --------------------------------------------------------------
% nfft = no_samples*2;
% % L
% hrtf_propL = squeeze(Obj.Data.IR(pos, 1, :));
% h_prop_smoothL = fft(hrtf_propL, nfft);
% h_prop_smoothL=h_prop_smoothL./h_prop_smoothL(f_1000hz); %normalizar em 1000Hz
% h_prop_smoothL = 20*log10(abs(h_prop_smoothL));
% % R
% hrtf_propR = squeeze(Obj.Data.IR(pos, 2, :));
% h_prop_smoothR = fft(hrtf_propR, nfft);
% h_prop_smoothR=h_prop_smoothR./h_prop_smoothR(f_1000hz); %normalizar em 1000Hz
% h_prop_smoothR = 20*log10(abs(h_prop_smoothR));
%
%
%
% %%%% PLOT %%%%
% figure()
% semilogx(freq(1:N/2), h_prop_smoothL(1:N/2), 'lineWidth', 1.5,  'Color', 'blue'); hold on
%
% semilogx(freq(1:N/2), h_prop_smoothR(1:N/2), '-.r', 'lineWidth', 1.5, 'Color', 'blue');
%
% legend('Esquerda', 'Direita', ...
%                'Location', 'best', 'lineWidth', 1.5)
%
% xlabel('Frequ�ncia [Hz]');
% ylabel('Amplitude [dB]');
%
% azim = num2str(out_pos(pos,1));
% elee = num2str(out_pos(pos,2));
% title(['DTF Individuo ', subject,', azimute ' azim '�, eleva��o ' elee '�.' ]);
%
% grid on
% axis tight
% set(gca,'FontSize',13)
%




% Davi R. Carvalho - Maio/2020
% Extract SOFA HRIRs to HeSuVi .wav file format
%% HeSuVi ANGLES
% HESUVI_TRACK_ORDER = ['FL-left', 'FL-right', 'SL-left', 'SL-right', 'BL-left', 'BL-right', 'FC-left', 'FR-right',
%                       'FR-left', 'SR-right', 'SR-left', 'BR-right', 'BR-left', 'FC-right']
% SPEAKER_ANGLES = {
%     'FL': 30,
%     'FR': -30, = 330
%     'FC': 0,
%     'BL': 150,
%     'BR': -150, = 210
%     'SL': 90,
%     'SR': -90 = 270
% }

%% INTERNAL FUNCTIONS
function HeSuViExport(Obj)
  Obj = sofaResample(Obj, 48000, 4096); %resample to 48kHz e aumentar tamanho do vetor de saida pra permitir low freq
  Obj = sofaAddLowFrequency(Obj);

  %% Find index of speaker angles
  azim = [30, 110, 150, 0, 330, 250, 210];
  elev = 0; n = 0;
  % Get index of measurements with the same directions
  Obj.SourcePosition = round(Obj.SourcePosition);
  for k = 1:length(azim)
      pos_idx = find(Obj.SourcePosition(:,1)==azim(k) & Obj.SourcePosition(:,2)==elev); %#ok<SAGROW>
      IR = Obj.Data.IR(pos_idx,:,:);
      for l = 1:size(IR, 2) % L/R channels
          n = n+1;
          y(:,n) = IR(:,l,:);
      end
  end
  %% Export to .wav
  y = y./max(abs(y(:)));
  y = y(:, [1,2,3,4,5,6,7,10,9,12,11,14,13,8]);
  filename = 'HRIR_140_48kHz.wav';
  Fs = Obj.Data.SamplingRate;
  audiowrite(filename,y,Fs)
end


function Obj = sofaResample(Obj, Fs, Nintp)
  N = ceil( (Fs/Obj.Data.SamplingRate) * size(Obj.Data.IR, 3) ); % length after resample
  if~Nintp
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
  %% Output
  Obj.Data.IR = IR;
  % update sampling rate
  Obj.Data.SamplingRate = Fs;
end

% Exemplo de reprodução de convolução em tempo real com filtros FIR
% DAVI ROCHA CARVALHO - janeiro/2021 @ UFSM
clear all; clc
addpath(genpath(pwd))
%% Carregar HRTF SOFA
Obj = SOFAload('individuo_141.sofa');
HRIRs = shiftdim(Obj.Data.IR, 2); % separar HRIRs 
positions = Obj.SourcePosition; % separar posições de fonte

azim = 40; % anti-horario
sourceDistance = 0.2;
elev = 0;


% Fazer parallax
samplingGrid = positions;
samplingGrid(:,2) = 90-samplingGrid(:,2);% Mod apenas para calculo
% sourceDistance = 0.3;
radius = 0.1; % raio da cabeça
[samplingGridParL, samplingGridParR] = supdeq_parallax(samplingGrid, sourceDistance, radius);





%% Carregar audio
buffer_size = 1024;
frameLength = 2^nextpow2(buffer_size);
local = 'Voice Sabine Short_edited.wav';
Audio = dsp.AudioFileReader(local, ...
                           'SamplesPerFrame', frameLength);

                       
%% Configura dispositivo de audio para reprodução
fs_Audio = Audio.SampleRate;      
deviceWriter = audioDeviceWriter('SampleRate', fs_Audio);
N_out = 2; % número de canais no audio de saída
setup(deviceWriter,zeros(Audio.SamplesPerFrame, N_out))


%% loop para reprodução em tempo real
azim = 40; % anti-horario
raio = 0.2;
elev = 0;
hazim = azim;
helev = elev;

% as duas linhas abaixo podem ser usadas para interromper a reprodução
release(Audio) 
release(deviceWriter)

% Inicializar output
Audio_out = zeros(buffer_size, N_out);

% Normalizar HRIR
HRIRs = HRIRs./max(abs(HRIRs(:)));
HRIRs = HRIRs * 0.09/raio;

% Inicializando filtros FIR
FIR_L = dsp.FIRFilter('NumeratorSource','Input port');
FIR_R = dsp.FIRFilter('NumeratorSource','Input port');

[posL, posR] = get_parallax_pos(positions, azim, elev, raio);    



%% Ambisonics ROOM
order = 7; % defina a ordem do decoder
idx_pos = get_pos(positions, azim, elev);
devices = positions(idx_pos, [1,2]).'; % axzimute e elevação dos "speakers"
dmtrx = audioexample.ambisonics.ambidecodemtrx(order, devices);
format = 'acn-sn3d';


% resample 
% desiredFs = fs_Audio;
[audio,fs] = audioread('0deg_066_Eigen_4th_Bformat_ACN_SN3D_44.wav');
% audio = resample(audio,desiredFs,fs);
% audiowrite('0deg_066_Eigen_4th_Bformat_ACN_SN3D_44.wav',audio,desiredFs);
audioDecoded = audioexample.ambisonics.ambidecode(audio, dmtrx, format);

%%
y = []; t=[];
release(Audio) 
release(deviceWriter)
fl = dsp.FIRFilter(HRIRs(:,posL,1).');
fr = dsp.FIRFilter(HRIRs(:,posR,2).');

while ~isDone(Audio)  
%     tic    
    % entrega um frame do audio a cada iteração
    audioIn = Audio(); 
    audioIn(:,2)=  audioIn(:,1);
    tic 
    Audio_out(:,1) = step(fl, audioIn(:,1));
    Audio_out(:,2) = step(fr, audioIn(:,2));
    t = [t toc];
    % Aplicar processamento de fato
%     tic
%     Audio_out(:,1) = FIR_L(audioIn(:,1), HRIRs(:,posL,1).'); % Esquerda
%     Audio_out(:,2) = FIR_R(audioIn(:,1), HRIRs(:,posR,2).'); % Direita
%     t = [t toc];
    
    % Reproduzir audio
    deviceWriter(Audio_out);    
    y = [y; Audio_out];
end
release(Audio) 
release(deviceWriter)
y = y./max(abs(y(:))); % normalizar saída
% filename = 'Sabine_Convo.wav';
% audiowrite(filename, y, fs_Audio)

mean(t)
%% conv ambisonics
nfft = length(audioDecoded) + length(y) - 1;
ambi = real(ifft(fft(audioDecoded, nfft) .* fft(y, nfft)));

sound(ambi, fs)

%% Funções internas -------------------------------------------------------
function[posL, posR] = get_parallax_pos(positions, azim, elev, sourceDistance)
% Fazer parallax
samplingGrid = positions;
samplingGrid(:,2) = 90-samplingGrid(:,2);% Mod apenas para calculo
% sourceDistance = 0.3;
radius = 0.1; % raio da cabeça
[samplingGridParL, samplingGridParR] = supdeq_parallax(samplingGrid, sourceDistance, radius);

% Pegar o indice da posição
posL = get_pos(samplingGridParL, azim, 90-elev);
posR = get_pos(samplingGridParR, azim, 90-elev);

% Pegar o indice da posição no grid original

% posL = get_pos(samplingGridParL, samplingGridParL(posL_tmp,1), samplingGridParL(posL_tmp,2));
% posR = get_pos(samplingGridParR, samplingGridParR(posR_tmp,1), samplingGridParR(posR_tmp,2));
end



function pos_idx = get_pos(positions, azim, elev)
    % Encontrar posição mais próxima à especificada
%     [~,pos_idx] = min(sqrt((positions(:,1) - azim).^2 + ...
%                            (positions(:,2) - elev).^2));
     pos_idx = dsearchn(positions(:, 1:2), [azim, elev]);

end

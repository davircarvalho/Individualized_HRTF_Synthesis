clear all; clc; tic
% DAVI ROCHA CARVALHO; ENG. ACUSTICA - UFSM; Fevereiro/2020
%% General INFO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ~INPUT SOFA FILE FROM THE CIPIC ARI AND ITA DATASETS,
% ~FIT TO THE CIPIC COORDINATE SYSTEM AND SAMPLE RATE
% ~PROCESS HRIR to CORRESPONDING DTF
%%% 
% ~Dados "remove_*.mat" vem da rotina "Anthropometry_Datasets.m"
%% PATHs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath([pwd, '\..\Functions']);
addpath([pwd, '\..\DADOS_TREINAMENTO']);
local = [pwd, '\..\Datasets\'];

% CIPIC
pathcipic = dir([local 'CIPIC\*.sofa']);
[~,idx_cipic] = natsortfiles({pathcipic.name});
pathcipic = pathcipic(idx_cipic, :);

% ARI
pathari = dir([local 'ARI\hrtf b_nh*.sofa']);
[~,idx_ari] = natsortfiles({pathari.name});
pathari = pathari(idx_ari, :);

% ITA
pathita = dir([local 'AACHEN\*.sofa']);
[~,idx_ita] = natsortfiles({pathita.name});
pathita =  pathita(idx_ita, :);

% 3D3A
path3d3a = dir([local '3D3A\Public-Data\Subject*\Subject*_HRIRs.sofa']);
[~,idx_3d3a] = natsortfiles({path3d3a.name});
path3d3a =  path3d3a(idx_3d3a, :);

% RIEC
pathriec = dir([local 'RIEC\*.sofa']);
[~,idx_riec] = natsortfiles({pathriec.name});
pathriec =  pathriec(idx_riec, :);

% TU Berlim 
pathtub_meas = dir([local 'HUTUBS\pp*_HRIRs_measured.sofa']);
[~,idx_tubmeas] = natsortfiles({pathtub_meas.name});
pathtub_meas = pathtub_meas(idx_tubmeas, :);


pathtub_sim = dir([local 'HUTUBS\pp*_HRIRs_simulated.sofa']);
[~,idx_tubsim] = natsortfiles({pathtub_sim.name});
pathtub_sim = pathtub_sim(idx_tubsim, :);


%% Options
% Defina quais datasets usar: {'cipic', 'ari', 'ita', '3d3a', 'riec', 'tub_meas', 'tub_sim'}, o
Datasets = {'ari','3d3a', 'ita', 'tub_sim'};
no_samples = 200; % Tamanho do vetor na saída (pós fft)
fs   = 44100;     % Taxa de amostragem 
fmin = 250;       % Frequencia min de corte para RI 
fmax = 18000;     % Frequencia max de corte para IR  
% Grid objetivo selecionado a partir do grid com menor número de posições
out_pos = select_best_grid(Datasets);
freq = linspace(0, fs-fs/no_samples, no_samples);

Datasets = {'ari', 'ita', 'tub_sim'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CIPIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(strcmp('cipic', Datasets))  
    disp('Processando CIPIC Dataset ...'); tic
    DTF_CIPIC = zeros( size(out_pos,1), 2, no_samples, length(pathcipic));
    for k = 1:length(pathcipic)
        CIPIC = SOFAload([pathcipic(k).folder '\' pathcipic(k).name]);
        % Process
        CIPIC_ok = process2unite(CIPIC, out_pos, fs, fmin, fmax);        
        DTF_CIPIC(:,:,:,k) = (abs(fft(CIPIC_ok.Data.IR, no_samples, 3)));
    end  
    DTF_CIPIC = permute(DTF_CIPIC, [3,4,1,2]);
    % remover dados sem antropometria e teste
    load('remove_CIPIC.mat');
    DTF_CIPIC(:,remove_CIPIC,:,:) = [];
    toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ARI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(strcmp('ari', Datasets))  %(hrtf b_nh*.sofa)
    disp('Processando ARI Dataset ...'); tic
    anthro_ari = load('anthro_ari.mat');
    no_subj_ari = length(pathari);    
    % Identificar todo o dataset (so o numero)
    for k = 1:no_subj_ari 
        id_hrtf_ari(k) = sscanf(pathari(k).name, 'hrtf b_nh%d');  
    end
    
    no_subj_ari = length(anthro_ari.id);
    DTF_ARI = zeros( size(out_pos,1), 2, no_samples, no_subj_ari);
    for k = 1:no_subj_ari
    % Encontrar indivíduos com antropometria
        idx = find((anthro_ari.id(k) - 3000) == id_hrtf_ari);
    % Carregar SOFA files
        ARI = SOFAload([pathari(idx).folder, '\',pathari(idx).name], 'nochecks');
    % Process
        ARI_ok = process2unite(ARI, out_pos, fs, fmin, fmax);
        DTF_ARI(:,:,:,k) = (abs(fft(ARI_ok.Data.IR, no_samples, 3)));
    end
    
    DTF_ARI = permute(DTF_ARI, [3,4,1,2]);
    % Remover HRTFs com antropometria INCOMPLETA
    load('remove_ARI.mat');
    DTF_ARI(:,remove_ARI,:,:) = [];
    toc
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ITA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(strcmp('ita', Datasets))
    disp('Processando ITA Dataset ...'); tic
    DTF_ITA = zeros( size(out_pos,1), 2, no_samples, length(pathita));
    for k = 1:length(pathita) 
        ITA = SOFAload([pathita(k).folder, '\',pathita(k).name], 'nochecks');
    
        %%% Transição de coordenadas cartesianas para esfericas
        for l = 1:length(ITA.SourcePosition)
            x = ITA.SourcePosition(l, 1);  
            y = ITA.SourcePosition(l, 2); 
            z = ITA.SourcePosition(l, 3);
            % new coordinates
            [az,elev,r] = cart2sph(x,y,z);
            azi=rad2deg(az); elev=rad2deg(elev);
            [azi,ele]   = nav2sph(azi,elev);
            % update coordinates
            ITA.SourcePosition(l, 1) = azi;
            ITA.SourcePosition(l, 2) = ele; 
            ITA.SourcePosition(l, 3) = round(r);
            % more metadata
            ITA.SourcePosition_Type = 'spherical';
            ITA.SourcePosition_Units = 'degree, degree, meter';              
        end       
                
        % Process
        ITA_ok = process2unite(ITA, out_pos, fs, fmin, fmax);
        DTF_ITA(:,:,:,k) = (abs(fft(ITA_ok.Data.IR, no_samples, 3)));        
    end
    % Remover indivíduos com dados inconsistentes
    DTF_ITA = permute(DTF_ITA, [3,4,1,2]);
    DTF_ITA(:,(14:15),:,:) = [];
    toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3D3A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(strcmp('3d3a', Datasets))  
    disp('Processando 3D3A Dataset ...'); tic
    % Encontrar HRIRs para devidas antropometria
    anthro_D3A = load('anthro_3D3A.mat');
    for k = 1:length(path3d3a)
        temp(k) = regexp(path3d3a(k).name, '\d+', 'match'); % get the id number 
        num(k) = str2num(temp{k});
    end
    
    % Carregar HRIRs
    DTF_D3A = zeros( size(out_pos,1), 2, no_samples, length(anthro_D3A.id));
    for k = 1:length(anthro_D3A.id)
        idx = find(num == anthro_D3A.id(k)); 
        D3A = SOFAload([path3d3a(idx).folder '\' path3d3a(idx).name], 'nochecks');       
        
        % Process
        D3A_ok = process2unite(D3A, out_pos, fs, fmin, fmax);
        DTF_D3A(:,:,:,k) = (abs(fft(D3A_ok.Data.IR, no_samples, 3)));       
    end  
    DTF_D3A = permute(DTF_D3A, [3,4,1,2]);  
    load('remove_D3A.mat');
    DTF_D3A(:,remove_D3A,:,:) = []; % remover por inconsistencias
    toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RIEC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(strcmp('riec', Datasets))  
    disp('Processando RIEC Dataset ...'); tic 
    DTF_RIEC = zeros( size(out_pos,1), 2, no_samples, length(pathriec));
    for k = 1: length(pathriec)
        RIEC = SOFAload([pathriec(k).folder '\' pathriec(k).name], 'nochecks');
        
        % Process
        RIEC_ok = process2unite(RIEC, out_pos, fs, fmin, fmax);
        DTF_RIEC(:,:,:,k) = (abs(fft(RIEC_ok.Data.IR, no_samples, 3)));             
    end  
    DTF_RIEC = permute(DTF_RIEC, [3,4,1,2]);    
    toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TU Berlim MEDIDO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(strcmp('tub_meas', Datasets))  
    disp('Processando TU Berlim (medido) Dataset ...'); tic
    DTF_TUBmeas = zeros( size(out_pos,1), 2, no_samples, length(pathtub_meas));
    for k = 1 : length(pathtub_meas)
        TUBmeas = SOFAload([pathtub_meas(k).folder '\' pathtub_meas(k).name], 'nochecks');                
        
        % Process
        TUBmeas_ok = process2unite(TUBmeas, out_pos, fs, fmin, fmax);
        DTF_TUBmeas(:,:,:,k) = (abs(fft(TUBmeas_ok.Data.IR, no_samples, 3)));  
    end  
    DTF_TUBmeas = permute(DTF_TUBmeas, [3,4,1,2]);
    load('remove_TUB.mat');
    DTF_TUBmeas(:,remove_TUB,:,:) = [];
    toc 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TU Berlim SIMULADO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(strcmp('tub_sim', Datasets))  
    disp('Processando TU Berlim (simulado) Dataset ...'); tic
    DTF_TUBsim = zeros( size(out_pos,1), 2, no_samples, length(pathtub_sim));
    for k = 1 : length(pathtub_sim)
        TUBsim = SOFAload([pathtub_sim(k).folder '\' pathtub_sim(k).name], 'nochecks');

        % Process
        TUBsim_ok = process2unite(TUBsim, out_pos, fs, fmin, fmax);
        DTF_TUBsim(:,:,:, k) = (abs(fft(TUBsim_ok.Data.IR, no_samples, 3)));
    end  
    DTF_TUBsim = permute(DTF_TUBsim, [3,4,1,2]);
    load('remove_TUB.mat');
    DTF_TUBsim(:,remove_TUB,:,:) = [];
    toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PFT Assembly %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DTF = [];% inicialização parcial
path_save = [pwd '\..\DADOS_TREINAMENTO\DTF'];
if any(strcmp('cipic', Datasets))
    DTF = cat(2, [DTF, DTF_CIPIC]);
    path_save = append(path_save, '_CIPIC');
end
if any(strcmp('ari', Datasets))
    DTF = cat(2, [DTF, DTF_ARI]);
    path_save = append(path_save, '_ARI');
end
if any(strcmp('ita', Datasets))
    DTF = cat(2, [DTF, DTF_ITA]);
    path_save = append(path_save, '_ITA');
end
if any(strcmp('3d3a', Datasets))
    DTF = cat(2, [DTF, DTF_D3A]);
    path_save = append(path_save, '_3D3A');
end
if any(strcmp('riec', Datasets))
    DTF = cat(2, [DTF, DTF_RIEC]);
    path_save = append(path_save, '_RIEC');
end
if any(strcmp('tub_meas', Datasets))
%     DTF = cat(2, [DTF, DTF_TUBmeas]);
%     path_save = append(path_save, '_TUBMEAS');
end
if any(strcmp('tub_sim', Datasets))
    DTF = cat(2, [DTF, DTF_TUBsim]);
    path_save = append(path_save, '_TUBSIM');
end
DTF = normDTF(DTF, freq);

save(path_save, 'DTF', 'out_pos', 'fs')
disp('Dados Salvos!')


%% plot
figure()
surf(DTF(:,:,500,2),'linestyle', 'none')


%% LOCAL FUCTIONS 
function Obj = process2unite(Obj, out_pos, fs, fmin, fmax)
    % Make same grid
    Obj = sofaFit2Grid(Obj, out_pos, 'adapt', 'Fs', fs);     
    % Normalize L/R balance and IR levels
    Obj = sofaNormalize(Obj);
    % filter
    Obj = sofaIRfilter(Obj, fmin, fmax);
    % HRTF -> DTF
    [Obj, ~] = SOFAhrtf2dtf(Obj);    
end


function DTF = normDTF(DTF, freq)
% f_idx = dsearchn(freq', 500); %posicao da sample em 500HZ
% DTF = DTF./DTF(f_idx,:,:,:);
DTF = 20*log10(DTF);
end
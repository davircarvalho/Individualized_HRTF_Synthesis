function out_pos = select_best_grid(Datasets)
%% Paths 
filepath = fileparts(which('select_best_grid.m'));
local = [filepath, '\..\Datasets\'];
% CIPIC
pathcipic = dir([local 'CIPIC\*.sofa']);

% ARI
pathari = dir([local 'ARI\hrtf b_nh*.sofa']);

% ITA
pathita = dir([local 'AACHEN\*.sofa']);

% 3D3A
path3d3a = dir([local '3D3A\Public-Data\Subject*\Subject*_HRIRs.sofa']);

% RIEC
pathriec = dir([local 'Banco RIEC\HRIRs\*.sofa']);

% TU Berlim Measured
pathtub_meas = dir([local 'HUTUBS\pp*_HRIRs_measured.sofa']);

% TU Berlim Simulated
pathtub_sim = dir([local 'HUTUBS\pp*_HRIRs_simulated.sofa']);

% Chedar
% pathchedar = dir('C:\Users\rdavi\Desktop\cheddar\sofacoustics.org\data\database\chedar\chedar_*_UV1m.sofa');

%% LOAD
if any(strcmp('cipic', Datasets))  
    CIPIC = SOFAload([pathcipic(1).folder '\' pathcipic(1).name], 'nochecks');
    POS{1} = CIPIC.SourcePosition;                         
end
if any(strcmp('ari', Datasets))  
    ARI = SOFAload([pathari(1).folder '\' pathari(1).name], 'nochecks');
    POS{2} = ARI.SourcePosition;                         
end
if any(strcmp('ita', Datasets))  
    ITA = SOFAload([pathita(1).folder '\' pathita(1).name], 'nochecks');
    POS{3} = ITA2sph(ITA.SourcePosition);                         
end
if any(strcmp('3d3a', Datasets))  
    D3A = SOFAload([path3d3a(10).folder '\' path3d3a(10).name], 'nochecks');
    POS{4} = D3A.SourcePosition;                         
end
if any(strcmp('riec', Datasets))  
    RIEC = SOFAload([pathriec(10).folder '\' pathriec(10).name], 'nochecks');
    POS{5} = RIEC.SourcePosition;                         
end
if any(strcmp('tub_meas', Datasets))  
    TUBmeas = SOFAload([pathtub_meas(1).folder '\' pathtub_meas(1).name], 'nochecks');
    POS{6} = TUBmeas.SourcePosition;                         
end
if any(strcmp('tub_sim', Datasets))  
    TUBsim = SOFAload([pathtub_sim(1).folder '\' pathtub_sim(1).name], 'nochecks');
    POS{7} = TUBsim.SourcePosition;                         
end
% if any(strcmp('chedar', Datasets))  
%     CHED = SOFAload([pathchedar(1).folder '\' pathchedar(1).name], 'nochecks');
%     POS{8} = CHED.SourcePosition;                         
% end

%% Compare
len = zeros(length(POS));
for k = 1:size(POS, 2)
    len(k) = length(POS{k}); %numero de posições em cada dataset
end
len(len==0) = inf; %evitar escolher datasets que não estão em uso 
[~,idx] = min(len);
out_pos = POS{idx};


%% Convert ITA to sph 
function pos = ITA2sph(ori_pos) 
    for l = 1:length(ori_pos)
        x = ori_pos(l, 1);  
        y = ori_pos(l, 2); 
        z = ori_pos(l, 3);
        % new coordinates
        [az,elev,r] = cart2sph(x,y,z);
        azi=rad2deg(az); elev=rad2deg(elev);
        [azi,ele]   = nav2sph(azi,elev);
        % update coordinates
        pos(l, 1) = azi;
        pos(l, 2) = ele; 
        pos(l, 3) = round(r);           
    end       
end
end
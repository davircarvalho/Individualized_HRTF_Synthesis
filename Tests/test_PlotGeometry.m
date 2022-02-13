% clear all; clc; tic
% DAVI ROCHA CARVALHO; ENG. ACUSTICA - UFSM; DEZEMBRO/2020
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

% VIKING
pathvik = dir([local 'VIKING\*.sofa']);
[~,idx_vik] = natsortfiles({pathvik.name}); % garantir que estão em ordem
pathvik = pathvik(idx_vik, :);


%% Output positions
res = 1;
ele=[-90:res:90 89:-res:-90 zeros(1,length(1:res:360-res))]';
azi=[zeros(length(-90:res:90),1); 180*ones(length(89:-res:-90),1); (1:res:360-res)'];
r = ones(length(azi),1);

out_pos = unique([azi, ele, r], 'rows');



%%
k = 10;
CIPIC = SOFAload([pathcipic(k).folder '\' pathcipic(k).name], 'nochecks');
ARI = SOFAload([pathari(k).folder, '\',pathari(k).name], 'nochecks');
ITA = ITA2spheric(SOFAload([pathita(k).folder, '\',pathita(k).name], 'nochecks'));
D3A = SOFAload([path3d3a(k).folder '\' path3d3a(k).name], 'nochecks');       
TUBmeas = SOFAload([pathtub_meas(k).folder '\' pathtub_meas(k).name], 'nochecks');                
TUBsim = SOFAload([pathtub_sim(k).folder '\' pathtub_sim(k).name], 'nochecks');
VIK = SOFAload([pathvik(k).folder '\' pathvik(k).name], 'nochecks');

%% Fix somethings

CIPIC = rmfield(CIPIC,'SourceView');
CIPIC = rmfield(CIPIC,'SourceUp');
TUBsim.ReceiverPosition = TUBsim.ReceiverPosition/10;
TUBmeas.ReceiverPosition = TUBmeas.ReceiverPosition/10;

% TUBsim.SourcePosition = inpt_pos;
% TUBsim.Data.IR = inpt_IR;
%%
close all; clc
plot_geometry(CIPIC,'CIPIC');
plot_geometry(ARI, 'ARI');
plot_geometry(ITA, 'AACHEN'); 
plot_geometry(D3A, '3D3A');
plot_geometry(TUBmeas, 'HUTUBS (measured)');
plot_geometry(TUBsim, 'HUTUBS (simulated)');







function plot_geometry(Obj, name)
    Obj.SourcePosition(:,3) = ones(size(Obj.SourcePosition(:,3),1),1);
    SOFAplotGeometry(Obj); hold on
    hFigure = gcf;
    view([35 20])
    
    ax = gca;
    axis(ax, 'off')
    title(name)
    legend('off')   
    set(gca, 'FontSize', 12)
    
    lm=1;
    xlim([-lm lm])
    ylim([-lm lm])
    zlim([-lm lm])
    
    filename = [pwd, '\Images\geometry_' name '.pdf' ];
    exportgraphics(hFigure,filename,'BackgroundColor','none','ContentType','vector')
end
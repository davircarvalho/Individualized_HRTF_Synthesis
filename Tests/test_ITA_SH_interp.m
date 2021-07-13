clear all; clc 
%%% Load data
local = [pwd '\..\Datasets\'];    
path = dir([local 'HUTUBS\pp*_HRIRs_simulated.sofa']);
Obj = SOFAload([path(1).folder, '\',path(1).name], 'nochecks');


%% ITA TOOLBOX
Obj_ita = SOFA2itaHRTF(Obj);
res = 1;
theta_ele=[-90:res:90 89:-res:-90 zeros(1,length(1:res:355))]';
phi_azi=[zeros(length(-90:res:90),1); 180*ones(length(89:-res:-90),1); (1:res:355)'];
r=ones(size(theta_ele));


pos = [theta_ele phi_azi, r];
% define interpolated positions
coords = itaCoordinates(size(pos,1));
coords.phi_deg = phi_azi;
coords.theta_deg = theta_ele+90;
coords.r = r;


%% SH INTER/EXTRAP
Obj_ita_extrp = Obj_ita.interp(coords, 'epsilon', 1e-8);
Obj_sofa_extrp = itaHRTF2SOFA(Obj_ita_extrp);



%%
plane = 'MagSagittal';
figure();
SOFAplotHRTF(Obj_sofa_extrp, plane); title('Interpolated')
axis tight

figure();
SOFAplotHRTF(Obj, plane); title('Reference')
axis tight


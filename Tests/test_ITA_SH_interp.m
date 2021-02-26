clear all; clc 
%%% Load data
local = [pwd '\..\Datasets\'];    
path = dir([local 'HUTUBS\pp*_HRIRs_simulated.sofa']);
Obj = SOFAload([path(1).folder, '\',path(1).name], 'nochecks');


%% ITA TOOLBOX
SOFAsave('temp.sofa', Obj);
Obj_ita = itaHRTF('sofa', 'temp.sofa');

res = 1;
theta=[-90:res:90 89:-res:-90 zeros(1,length(1:res:355))]';
phi=[zeros(length(-90:res:90),1); 180*ones(length(89:-res:-90),1); (1:res:355)'];
r=ones(size(theta));

% define interpolated positions
coords = itaCoordinates([r, r, r], 'cart');
coords.phi_deg = phi;
coords.theta_deg = theta;
coords.r = r;


%% SH INTER/EXTRAP
Obj_ita_extrp = Obj_ita.interp(coords);



%% save SOFA 
Obj_sofa_extrp = SOFAgetConventions('SimpleFreeFieldHRIR');

% IR
Obj_sofa_extrp.Data.SamplingRate = Obj_ita_extrp.samplingRate;
irL = Obj_ita_extrp.getEar('L');
irR = Obj_ita_extrp.getEar('R');
Obj_sofa_extrp.Data.IR = zeros(size(irL.time, 2), 2, size(irL.time, 1));
Obj_sofa_extrp.Data.IR(:,1,:) = (irL.time).';
Obj_sofa_extrp.Data.IR(:,2,:) = (irR.time).';

% Coordinates
Obj_sofa_extrp = ita_itaHRTFCoordinates2SOFACoordinates(Obj_ita_extrp, Obj_sofa_extrp);

% Meta
filename = (['ita.sofa']);
Obj_sofa_extrp = SOFAupdateDimensions(Obj_sofa_extrp);
SOFAsave(filename, Obj_sofa_extrp);
disp('Done!')

%%
plane = 'MagSagittal';
figure();
SOFAplotHRTF(Obj_sofa_extrp, plane); title('Interpolated')
axis tight

figure();
SOFAplotHRTF(Obj, plane); title('Reference')
axis tight


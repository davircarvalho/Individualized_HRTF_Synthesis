clear all; clc

addpath('../Functions')
Obj = SOFAload('../Datasets/CIPIC/Subject_003.sofa');

result = SOFA2itaHRTF(Obj);
Obj_out = itaHRTF2SOFA(result);


%% Plot 
plane1 = 'maghorizontal';
figure()
SOFAplotHRTF(Obj,plane1); title(['Reference - ' plane1]);
axis tight
xlim([0 2e4])

figure()
SOFAplotHRTF(Obj_out,plane1); title(['Output - ' plane1]);
axis tight
xlim([0 2e4])


%% 
%%
pos = Obj.SourcePosition; 
coordinates = itaCoordinates(size(pos,1));
coordinates.phi_deg = pos(:,1);
coordinates.theta_deg = pos(:,2)+90;
coordinates.r = pos(:,3);   
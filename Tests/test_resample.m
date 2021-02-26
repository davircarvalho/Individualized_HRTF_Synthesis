clear all; clc
addpath([pwd, '\..\Functions']);
local = [pwd, '\..\Datasets\3D3A\Public-Data\Subject17\Subject17_HRIRs.sofa'];
Obj = SOFAload(local);


%%
Obj.Data.IR = Obj.Data.IR ./ max(abs(Obj.Data.IR(:)));
N1 = size(Obj.Data.IR, 3);
Obj2 = sofaResample(Obj, 44100);
N2 = size(Obj2.Data.IR, 3);

fs1 = Obj.Data.SamplingRate;
fs2 = Obj2.Data.SamplingRate;

f1 = linspace(0, fs1-fs1/N1, N1);
f2 = linspace(0, fs2-fs2/N2, N2);
IR1 = shiftdim(Obj.Data.IR, 2);
IR2 = shiftdim(Obj2.Data.IR, 2);


%% plot
subplot(211)
tx1 = (0:N1-1)/fs1; 
tx2 = (0:N2-1)/fs2; 
ir1 = IR1(:,1,1) ./ max(abs(IR1(:,1,1)));
ir2 = IR2(:,1,1) ./ max(abs(IR2(:,1,1)));

plot(tx1, ir1); hold on 
plot(tx2, ir2); hold on 
axis tight
legend('original', 'resample', 'location', 'best')
xlabel('tempo (s)')
subplot(212)
semilogx(f1, db(abs(fft(ir1, N1)))); hold on
semilogx(f2, db(abs(fft(ir2, N2))));
xlim([0, 2e4])
legend('original', 'resample', 'location', 'best')
xlabel('frequencia (Hz)')

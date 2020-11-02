%%
clc
addpath([pwd, '\Functions']);
local = [pwd, '\Datasets\3D3A\Public-Data\Subject17_HRIRs.sofa'];
Obj = SOFAload(local);


%%
Obj.Data.IR = Obj.Data.IR ./ max(abs(Obj.Data.IR(:) ));
clc
N = size(Obj.Data.IR, 3);
fs1 = Obj.Data.SamplingRate;
Obj2 = sofaResample(Obj, 44100, N);

fs2 = Obj2.Data.SamplingRate;

f1 = linspace(0, fs1-fs1/N, N);
f2 = linspace(0, fs2-fs2/N, N);
IR1 = shiftdim(Obj.Data.IR, 2);
IR2 = shiftdim(Obj2.Data.IR, 2);


%% plot
subplot(211)
tx1 = 0:1/fs1:N/fs1; tx1(end) = [];
tx2 = 0:1/fs2:N/fs2; tx2(end) = [];
plot(tx1, IR1(:,1,1)); hold on 
plot(tx2, IR2(:,1,1)); hold on 
axis tight
legend('original', 'resample', 'location', 'best')

subplot(212)
semilogx(f1, db(abs(fft(IR1(:,1,1), N)))); hold on
semilogx(f2, db(abs(fft(IR2(:,1,1), N))));
xlim([0, 2e4])
legend('original', 'resample', 'location', 'best')
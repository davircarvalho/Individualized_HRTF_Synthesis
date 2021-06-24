% clear all; clc
addpath([pwd, '\..\Functions']);
% local = [pwd, '\..\Datasets\3D3A\Public-Data\Subject17\Subject17_HRIRs.sofa'];
% local = [pwd, '\..\Datasets\ARI\hrtf b_nh8.sofa'];
local = [pwd, '\..\Datasets\CIPIC\subject_008.sofa'];

Obj = SOFAload(local);

%% Filter
fmin = 250;
fmax = 18000;
Obj2 = sofaIRfilter(Obj, fmin, fmax);

N = size(Obj.Data.IR, 3);
fs1 = Obj.Data.SamplingRate;
f = linspace(0, fs1-fs1/N, N);
IR1 = shiftdim(Obj.Data.IR, 2);
IR2 = shiftdim(Obj2.Data.IR, 2);

%% PLOT -------------------------------------------------------------------
figure;
% Time
subplot(311)
tx = (0:N-1)/fs1; 
ir1 = IR1(:,1,1);       
ir2 = IR2(:,1,1);
plot(tx, ir1); hold on 
plot(tx, ir2); hold on 
axis tight
legend('Original', 'filtered', 'location', 'best')
xlabel('tempo (s)')

% Frequency
subplot(312)
semilogx(f, db(abs(fft(ir1, N)))); hold on
semilogx(f, db(abs(fft(ir2, N))));
xlim([0, 2e4])
legend('original', 'filtered', 'location', 'best')
xlabel('frequencia (Hz)')

% Phase
subplot(313)
semilogx(f, imag(fft(ir1, N))); hold on
semilogx(f, imag(fft(ir2, N)));
xlim([0, 2e4])
legend('original', 'filtered', 'location', 'best')
xlabel('frequencia (Hz)')
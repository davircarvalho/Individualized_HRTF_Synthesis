clear all; clc
addpath([pwd, '\..\Functions']);
local = [pwd, '\..\Datasets\3D3A\Public-Data\Subject17\Subject17_HRIRs.sofa'];
% local = [pwd, '\..\Datasets\ARI\hrtf b_nh8.sofa'];
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


% plot
figure
subplot(211)
tx1 = (0:N1-1)/fs1; 
tx2 = (0:N2-1)/fs2; 
ir1 = IR1(:,1,1);       
ir2 = IR2(:,1,1);

plot(tx1, ir1); hold on 
plot(tx2, ir2); hold on 
axis tight
legend('Original', 'Resample', 'location', 'best')
xlabel('tempo (s)')
subplot(212)
semilogx(f1, db(abs(fft(ir1, N1)))); hold on
semilogx(f2, db(abs(fft(ir2, N2))));
xlim([0, 2e4])
legend('original', 'resample', 'location', 'best')
xlabel('frequencia (Hz)')






%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs1 = 96000;
t1 = 0:1/fs1:1;
x = chirp(t1, 0, 1, fs1/2, 'quadratic');
figure
spectrogram(x,kaiser(256,15),220,412,fs1,'yaxis')

%%
fs2 = 44100;
[p,q] = rat(fs2/fs1);
y = resample(x,p,q);
spectrogram(y,kaiser(256,15),220,412,fs2,'yaxis')

%%
fs2 = 44100;
t2 = 0:1/fs2:1;

y = resample_ya(x, fs1,fs2);
figure
spectrogram(y,kaiser(256,15),220,412,fs2,'yaxis')



function IR = resample_ya(x, Fs_in, Fs_out, Nintp)
N = ceil((Fs_out/Fs_in) * size(x, 3)); % length after resample
if nargin<4 || Nintp<N
    Nintp = N;
end
zpad = zeros((Nintp - N), 1);

N_sig = length(x);
%% options
tx = (0:N_sig-1)/Fs_in;
[p,q] = rat(Fs_out / Fs_in, 0.0001);
normFc = .98 / max(p,q);
order = 512 * max(p,q);
beta = 12;
%%% Cria um filtro via Least-square linear-phase FIR filter design
lpFilt = firls(order, [0 normFc normFc 1],[1 1 0 0]);
lpFilt = lpFilt .* kaiser(order+1,beta)';
lpFilt = lpFilt / sum(lpFilt);
% multiply by p
lpFilt = p * lpFilt;
% Actual Resample

IRpre = resample(x,p,q,lpFilt);
% IRpre = resample(x, tx, Fs_out, p,q, 'spline');
IR = [IRpre; zpad];
%% Output
% norm = max(abs(Obj.Data.IR));
IR = IR.* q/p;
end


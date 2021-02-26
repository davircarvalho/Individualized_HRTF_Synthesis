clear all; clc; close all;
Obj = SOFAload([pwd, '\..\Datasets\CIPIC\subject_003.sofa']);
addpath(genpath([pwd, '\..\Functions']));


%% preprocess
fs = Obj.Data.SamplingRate;
IR = shiftdim(Obj.Data.IR, 2);
N = size(IR, 1);
freq_vec = (0:N/2-1)*fs/N;


%%% prepare extension
fmin = 15; % minimum freq
N_ext = ceil(fs/fmin); % minimum length to have the fmin
if N_ext <= N
    N_ext = N;
    freq_vec_ext = freq_vec;
else
    freq_vec_ext = (0:N_ext/2-1)*fs/N_ext;
end

f500Hz = dsearchn(freq_vec_ext.', 500); % idx at 500Hz


%% interp 
% figure()
ir_interp = zeros(size(IR, 2), size(IR, 3), N_ext);
for k = 1:size(IR, 2)
    for l = 1:size(IR, 3)
        time = [IR(:,k,l); zeros(ceil(N_ext - N), 1)];
        mag = fft(time);
        mag_interp = mag;
        
        % interp baixas freqs
        x = [freq_vec_ext(2),    freq_vec_ext(f500Hz:f500Hz+1)];
        xq = freq_vec_ext(2:f500Hz);

        y_mag = [mag(f500Hz); mag(f500Hz:f500Hz+1)];
        
        mag_interp(2:f500Hz) = interp1(x, y_mag, xq, 'makima');

        % back to time domain
        ir_interp(k,l,:) = real(ifft(mag_interp, N_ext, 'symmetric'));
        

        %plot
%         semilogx(freq_vec_ext, db(abs(mag(1:N_ext/2)))); hold on
%         semilogx(freq_vec_ext, db(abs(mag_interp(1:N_ext/2)))); hold on 
%         hghg = fft(squeeze(ir_interp(k,l,:)));
%         semilogx(freq_vec_ext, db(abs(hghg(1:N_ext/2)))); hold off
%         legend('original', 'interp', 'ir_fft')
%         pause(1)
    end
end


%% Normalize
ir_interp = ir_interp./max(abs(ir_interp(:))) .* max(abs(IR(:)));

% OUTPUT
Obj_out = Obj;
Obj_out.Data.IR = ir_interp;


%% Plot time                       
% figure()
% tx = 0:1/fs:(N-1)/fs;
% tx_ext = 0:1/fs:(N_ext-1)/fs;
% 
% ch = 2;
% k=20;
% 
% plot(tx, IR(:,k,ch)); hold on
% plot(tx_ext(1:N), squeeze(ir_interp(k,ch,1:N))); hold off
% plot(tx_ext, squeeze(ir_interp(k,ch,:))); hold off
% 
% legend('original', 'interp')
% xlabel('Time (ms)')
% ylabel('Amplitude')
%             
% 
% %%%
% ch = 1;
% k = 1000;
% intr = db(abs(fft(squeeze(ir_interp(k,ch,:)))));
% ori = db(abs(fft(squeeze(Obj.Data.IR(k,ch,:)))));
% semilogx(freq_vec, ori(1:N/2)); hold on
% semilogx(freq_vec_ext, intr(1:N_ext/2));
% legend('original', 'interp')
% xlabel('Time (ms)')
% ylabel('Amplitude')
% 
% 
% %%% Plot ITD
% itd = sofaGetITD(Obj);
% itd2= sofaGetITD(Obj_out);
% figure()
% plot(itd); hold on; plot(itd2);
% legend('original', 'interp')

 
 
%% 
t = 0:1/1e3:10;
fo = 0;
f1 = 500;
y = chirp(t,fo,t(end),f1,'linear',0,'complex');

semilogx(angle(fft(y)))








 
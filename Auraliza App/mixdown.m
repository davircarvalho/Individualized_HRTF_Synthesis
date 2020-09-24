%% load
[x.f1, fs.f1] = audioread('GothicChurch.wav');
[x.f2,  fs.f2] = audioread('motetus.wav');
[x.f3, fs.f3] = audioread('tenor.wav');
[x.f4, fs.f4] = audioread('triplum.wav');

audF = fft(x.f4);
f1f = fft(x.f1, length(x.f4));
y = f1f.*audF;

%% avisa se fs diferentes
if fs.f1 ~= fs.f2 || ...
   fs.f2 ~= fs.f3 || ...
   fs.f3 ~= fs.f4
   warning('Taxas de amostragem diferentes, vai dar merda!')
end
%% process
len(1) = length(x.f1);
len(2) = length(x.f2);
len(3) = length(x.f3);
len(4) = length(x.f4);

[~, idx_max] = max(len);

%%
for k = 1:length(len)
    zpad = zeros(abs(len(k) - len(idx_max)), size(x.(['f' num2str(k)]), 2)); 
    x.(['paded_f' num2str(k)]) = [x.(['f' num2str(k)]); zpad];
end

%% Concatena
% faz separado para caso algum dos inputs não seja mono
y = [];
for k = 1:length(len)
    y = cat(2, y, x.(['paded_f' num2str(k)]));
end

%% Export

audiowrite('coral4ch.wav', y, fs.f1)
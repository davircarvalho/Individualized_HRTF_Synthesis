function sd = spec_dist(msrd,ref,fs,fmin,fmax)
% This function calculates the Spectral Distortion between a reference
% 'ref' and a measured data 'msrd' with length N. 
% Both vectores must be in TIME domain
% fmin and fmax limitate the frequency range to be compared
%% Equals the size of both vectors
M = length(msrd);           % stores in 'N' the length of msrd
N = length(ref);           % stores in 'M' the length of ref
% equals the size of the vectores in case that they are different
if M > N
    msrd(M) = 0;
    N = M;
end
if N > M
    ref(N) = 0;
end
%% Create frequency vector to find value in 1000 Hz
f = linspace(0,fs*(N-1)/N,N);
f_1000hz = dsearchn(f',1000);
fmin = dsearchn(f',fmin);
fmax = dsearchn(f',fmax);
%% Normalizes vectors in frequency domain
msrd = 2*fft(msrd, N); msrd = msrd/msrd(f_1000hz);
ref = 2*fft(ref, N);   ref = ref/ref(f_1000hz);

%% Calculates the Spectral Distortion
sd = sqrt((1/(fmax-fmin+1))*sum((20*log10(abs(msrd(fmin:fmax)./ref(fmin:fmax)))).^2));
end
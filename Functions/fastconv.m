function y = fastconv(x,h,trunc)
h = h(:);
x = x(:);

nfft = length(x)+length(h)-1;
y = ifft(fft(x,nfft).*fft(h,nfft));

if exist('trunc')
    y = y(1:length(x));
end
y = y(:);
end
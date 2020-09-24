function [Obj, itd] = sofaMinimunPhaseUP(Obj)
% input: objeto sofa
% output: objeto sofa com hrirs em fase minima (ITD é salvo em Delay)

% To find the minimum phase HRIR using Hilbert transform
% The minimum phase response and magnitude spectrum are connected
% via the Hilbert transform
% phaseRespMin = imaginary(HilbertTransform(-ln(magnitude of hSignal)))
hSignal = fft(shiftdim(Obj.Data.IR, 2));

hMag = abs(hSignal);            
hLnMag = log(hMag);            
hHilbert = hilbert(-hLnMag);            
hMinPhase = imag(hHilbert);

hMin = ifft((hMag.*exp(1i*hMinPhase)));
hMin = real(hMin);

itd_x = sofaGetITD(Obj);

offset = 30;
for k = 1:length(itd_x)
    pos = Obj.SourcePosition(k,:);
    if pos(1) >= 180 || (pos(2) < 0 && pos(1) >= 180)
        itd(k, 1) = round(itd_x(k) + offset);
        itd(k, 2) = offset; 
    else
        itd(k, 1) = offset;
        itd(k, 2) = round(itd_x(k) + offset);
    end
end

%% Output
% IR
Obj.Data.IR = shiftdim(hMin, 1);
end
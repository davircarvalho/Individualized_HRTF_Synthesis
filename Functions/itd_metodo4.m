% INPUT: respostas impulsivas de ambas as orelhas e vetor de frequencias
% Calculo do ITD por meio da relação de fase entre as HRIRs (IPD)
function [itd, naz, nel] = itd_metodo4(hrir_l, hrir_r, freqs)
if nargin < 3
    fs = 44100;
    N = length(hrir_l);  
    f = linspace(0, fs-fs/N, N);
    freqs = f(1:N/2);
end

% Azimutes
naz = [-80, -65, -55 -45:5:45, 55, 65, 80]; 
% Elevações
nel = -45 + 5.625*(0:49); 


end
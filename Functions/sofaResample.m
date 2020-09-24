function Obj = sofaResample(Obj, Fs, Nintp)
% Muda resolução aparente de objeto SOFA para valor especificado
% e faz zero padding para 2^nextpow2(N)

% Davi R. Carvalho @UFSM - Engenharia Acustica - Setembro/2020

%   Input Parameters:
%    Obj:        Objeto de HRTFs SOFA a ser modificado 
%    Fs:         Taxa de amostragem Objetivo
%    Nint (opcional):       Comprimento do vetor na saida, (adiciona zeros no final)

%   Output Parameters:
%     Obj_out:   Objeto de HRTFs SOFA com a taxa de amostragem Fs
%
% Matlab 2020a
%% Resample
N = ceil( (Fs/Obj.Data.SamplingRate) * size(Obj.Data.IR, 3) ); % length after resample
if nargin == 3 && Nintp < N
    Nintp = 2^nextpow2(N); % output length
else
    Nintp = N; % doesn't zero pad 
end
zpad = zeros((Nintp - N), 1);
% options
[p,q] = rat(Fs / Obj.Data.SamplingRate);
normFc = .98 / max(p,q);
order = 256 * max(p,q);
beta = 12;
%%% Cria um filtro via Least-square linear-phase FIR filter design
lpFilt = firls(order, [0 normFc normFc 1],[1 1 0 0]);
lpFilt = lpFilt .* kaiser(order+1,beta)';
lpFilt = lpFilt / sum(lpFilt);
% multiply by p
lpFilt = p * lpFilt;
% Actual Resample
for k = 1:size(Obj.Data.IR, 1)
    for l = 1:size(Obj.Data.IR, 2)
        IRpre(k, l, :) = resample(Obj.Data.IR(k, l, :),p,q,lpFilt);
        IR(k, l, :) = [squeeze(IRpre(k, l, :)); zpad];
    end 
end
%% Output
Obj.Data.IR = IR;
% update sampling rate
Obj.Data.SamplingRate = Fs;
end
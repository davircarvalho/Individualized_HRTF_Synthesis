function [hminL, hminR] = minphase(hL, hR, itd, azi)
% Calcula a fase mínima a partir da magnitude apenas
% Aplica a fase de excesso (ITD) sobre a RI [ITA TOOLBOX]

%%% INPUT %%%%
% hL, hR: magnitudes para orelha esquerda e direita (não espelhada e linear)
% itd: diferença iteraural de tempo (opcional)
% azi: azimute referente ao determinado itd (opcional)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pL = itaAudio; pR = pL;
pL.freq = hL;
pR.freq = hR;

%% Minimal phase
pminL = ita_minimumphase(pL, 'window', '0',...
                               'cutoff', 'false');                     
pminR = ita_minimumphase(pR, 'window', '0',...
                               'cutoff', 'false');

%% All pass phase (ITD)
if nargin <= 2
    itd = 0;
    azi = 1;
end
offset = 0;
if azi >= 0 
    pminL = ita_time_shift(pminL, itd+offset, 'samples');
    pminR = ita_time_shift(pminR, offset, 'samples');
else
    pminL = ita_time_shift(pminL, offset, 'samples');
    pminR = ita_time_shift(pminR, itd+offset, 'samples');
end

%% verificando se é fase minima
% isminphase(pminL.time)

%% save vector 
hminL = (pminL.time);
hminR = (pminR.time);

end
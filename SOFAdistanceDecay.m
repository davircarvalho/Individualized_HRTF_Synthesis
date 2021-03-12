function Obj = SOFAdistanceDecay(Obj, pos, varargin)
% Davi Rocha Carvalho MARCO/2021 @ UFSM
% Inputs:
%     pos = [azim, elev, radius] --> sound source positions
%     Obj --> SOFA SimpleFreeField HRTF
%
% Opcional:
%     'DistDecay'
%     'AirAbsorption'


%% Parse arguments
defaultAirAbs = 'AirAbsorption';
acceptedAirAbs = {'AirAbsorption', 'noAirAbsorption'};
checkAirAbs = @(x) any(validatestring(x, acceptedAirAbs));

defautDecay = 'DistDecay';
acceptedDecay = {'DistDecay', 'noDistDecay'};
checkDecay = @(x) any(validatestring(x, acceptedDecay));


p = inputParser;
addRequired(p,'Obj', @isstruct);
addRequired(p,'pos', @ismatrix);
addOptional(p,'decay', defautDecay,checkDecay)
addOptional(p,'airAbs', defaultAirAbs,checkAirAbs)
parse(p, Obj, pos, varargin{:})


%% General
ir = shiftdim(Obj.Data.IR, 2);
N = size(Obj.Data.IR, 3);
fs = Obj.Data.SamplingRate;
freq = linspace(0, fs-fs/N, N);
[~, Obj] = SOFAgetITD(Obj, 'samples');   
L_ear_pos = Obj.ReceiverPosition(1,:);
R_ear_pos = Obj.ReceiverPosition(2,:);


%% Calcular absorcao do ar por frequencia
if strcmp(p.Results.airAbs, defaultAirAbs)
    air_abs = zeros(length(freq), 1);
    for k = 1:length(freq)
        [~, air_abs(k,1), ~, ~] = air_absorption(freq(k));                    
    end
    air_abs(:,2) = air_abs(:,1);  
end


%% Processamento
for k = 1:size(pos, 1)
    %% Aplicar absorcao do ar em funcao da distancia
    if strcmp(p.Results.airAbs, defaultAirAbs)
        idx_pos = dsearchn(Obj.SourcePosition(:,1:2), [pos(k,1), pos(k,2)]);
        IR = [ir(:,idx_pos,1),...
              ir(:,idx_pos,2)];
        LR = 20*log10(abs(fft(IR))) - air_abs*pos(k,3);
        LR = real(ifft(get_min_phase(LR(1:N/2, :),'log','nonsymmetric')));
        % Devolver itd removido na fase minima
        ir(:,idx_pos,1) = circshift(LR(:,1), Obj.Data.Delay(idx_pos, 1), 1);
        ir(:,idx_pos,2) = circshift(LR(:,2), Obj.Data.Delay(idx_pos, 2), 1);
         % plot(freq, air_abs*pos(k,3)); hold on
    end
        %% Calcular decaimento com a distancia
    if strcmp(p.Results.decay, defautDecay)
        % Distancia da fonte para cada orelha 
        [tx, ty, tz] = sph2cart(pos(k,1), pos(k,2), pos(k,3));       
        Ldist = sqrt((L_ear_pos(1) - tx)^2 + (L_ear_pos(2) - ty)^2 + (L_ear_pos(3) - tz)^2);
        Rdist = sqrt((R_ear_pos(1) - tx)^2 + (R_ear_pos(2) - ty)^2 + (R_ear_pos(3) - tz)^2);
        % Fator de normalizacao
        DistNorm = [1/Ldist, 1/Rdist];
        ir(:,idx_pos,1) = ir(:,idx_pos,1).*DistNorm(1);
        ir(:,idx_pos,2) = ir(:,idx_pos,2).*DistNorm(2);
    end
end
Obj.Data.IR = shiftdim(ir,1);   
end

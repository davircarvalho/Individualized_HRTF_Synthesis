function sd = sofaSpecDist(Obj_sim, Obj_ref, fmin, fmax, varargin)
% This function calculates the Spectral Distortion between a reference
% 'ref' and a measured data 'msrd' with length N for SOFA object. 
% Both vectores must be in TIME domain
% fmin and fmax limitate the frequency range to be compared

% OPTIONAL ARGUMENTS:
%                    'freq':      will result the mean of the SD for all
%                                 positions at the first collumn and the std
%                                 deviation on the second collumn

%      (default)     'posi':      results the SD for each position,
%                                 the rms value per frequency vector is taken.

%% Parse Arguments
    % Método de processamento
    defaultMethod = 'posi';
    validMethods = {'posi','freq'};
    checkMethod = @(x) any(validatestring(x,validMethods));

    %Verificar entradas
    p = inputParser;
    addRequired(p,'Obj_sim', @isstruct);
    addRequired(p,'Obj_ref', @isstruct);
    addRequired(p,'fmin', @isnumeric);
    addRequired(p,'fmax', @isnumeric);
    addOptional(p,'method', defaultMethod, checkMethod);

    parse(p, Obj_sim,Obj_ref,fmin,fmax,varargin{:})


%% Distorção espectral ----------------------------------------------------
    fs = Obj_ref.Data.SamplingRate;
    % chama calculo por posicoes
    if length(Obj_ref.SourcePosition) ~= length(Obj_sim.SourcePosition)
        warning(['Mismatch -- Number of source positions is different! ', ...
                  'It will be evaluated only the positions corresponding to the '...
                  'coordinates in the object', ...
                  'with the least number of measurements. '...
                  'It will be applied nearest interpolation method to the ',...
                  'other SOFA file.'])
        
       if length(Obj_ref.SourcePosition) > length(Obj_sim.SourcePosition)
           Obj_ref = sofaFit2Grid(Obj_ref, Obj_sim.SourcePosition, 'adapt');
       else
           Obj_sim = sofaFit2Grid(Obj_sim, Obj_ref.SourcePosition, 'adapt');
       end
    end
    for k = 1:size(Obj_ref.Data.IR,1) % numero de posições 
        sim = squeeze(Obj_sim.Data.IR(k, 2, :));
        ref = squeeze(Obj_ref.Data.IR(k, 2, :));
        sd(:,k) = spec_dist(sim, ref, fs, fmin, fmax, p.Results.method);   
    end

%% Caso 'freq' mode -> mean and std deviation
    if strcmp('freq', p.Results.method)
        sd_mean = mean(sd, 2); % media de todas as posicoes
%         sd_std = std(sd,0,2); % desvio padrão
        clear sd
%         sd = [sd_mean, sd_std];
        sd = sd_mean;
    end
end




%%%%%%%%%%%%%%%%%%% INTERNAL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sd = spec_dist(msrd,ref,fs,fmin,fmax, method)
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
    f = linspace(0, fs-fs/N, N);
    f_500hz = dsearchn(f',500);
    fmin = dsearchn(f',fmin);
    fmax = dsearchn(f',fmax);
%% Normalizes vectors in frequency domain
    msrd = fft(msrd, N);
    msrd = msrd./msrd(f_500hz); % normalize @ 500Hz
    ref  = fft(ref, N);
    ref = ref./ref(f_500hz);    % normalize @ 500Hz

%% Calculates the Spectral Distortion
    switch method
        case 'posi'
            % por posição
            sd = sqrt((1/(fmax-fmin+1))*sum((20*log10(abs(msrd(fmin:fmax))./abs(ref(fmin:fmax)))).^2));
        case 'freq'
            % por frequência
            sd = 20*log10(abs(msrd(fmin:fmax)./ref(fmin:fmax)));
    end
end
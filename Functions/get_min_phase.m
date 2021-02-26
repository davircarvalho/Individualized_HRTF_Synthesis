function Hmin = get_min_phase(H, varargin)
% Cálculo de fase mínima a partir da magnitude apenas. get_min_phase opera 
% por colunas, para H(M x N) em que N é o número respostas de magnitude.
%
% Entradas: 
%           H: magnitude
% 
% Flags: 
%           'linear':        especifica se a magnitude de entrada é linear (padrão)
%           'log':           especifica se a magnitude de entrada é logaritimica
%
%           'symmetric':     especifica que a magnitude é simétrica
%           'nonsymmetric':  especifica que a magnitude é simétrica (padrão)
% 
% Saída:    Hmin:            espectro complexo de fase mínima, calculado a partir da
%                            Transformada de Hilbert.
%       
% Uso:                       Hmin = minimum_phase(H, 'linear', 'simetrica')

%% Parse arguments
defautScale = 'linear';
acceptedScales = {'linear', 'log'};
checkScale = @(x) any(validatestring(x, acceptedScales));

defaultSymmetry = 'symmetric';
acceptedSymmetries = {'symmetric', 'nonsymmetric'};
checkSymmetry = @(x) any(validatestring(x, acceptedSymmetries));

p = inputParser;
addRequired(p,'H',@ismatrix);
addOptional(p,'scale',defautScale,checkScale)
addOptional(p,'symmetry',defaultSymmetry,checkSymmetry)
parse(p, H, varargin{:})

%% Preprocess
if strcmp(p.Results.scale, 'log')
    H = 10.^(H./20); % back to linear
end

if strcmp(p.Results.symmetry, 'nonsymmetric')
    H = [H; flip(H)]; % back to double sided spectrum
end

%% Get minimum_phase
phi_min = imag(hilbert(-(log(abs(H)+eps)))); % eps makes avoids log(0) = -inf
% Filtro inverso complexo (cria fase)
Hmin = abs(H).*exp(1i*(phi_min));
end
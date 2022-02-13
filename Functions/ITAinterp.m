function cThis = ITAinterp(this,varargin)
% function this = interp(varargin)
%
% Function to calculate HRTFs for arbitrary field points using a N-th order
% spherical harmonics (SH) interpolation / range extrapolation, as described in [1],
% SH expansion coefficients are calculated by means of a least-squares
% approach with Tikhonov regularization
%
% Function may also be used for spatial smoothing of HRTF using
% the method described in [2]. As field input use the original
% measurement grid and set the desired order of the SH matrix /
% truncation order.
%
% INPUT:
%     varargin{1}      ...  itaCoordinates object (required)
%                           varargin{1}.phi: desired azimuth angles for HRTF interpolation [0 2*pi)
%                           varargin{1}.theta: desired zenith angles for HRTF interpolation [0 pi]
%                           varargin{1}.r: (optional) desired radius used for range extrapolation in [m],
%                                    set to 1 if no range extrapolation is required
%     order            ...  order of spherical harmonics matrix (default: 50)
%     epsilon          ...  regularization coefficient (default: 1e-8)
%     shiftToEar            to a shift to approximate ear position to
%                           improve sh transformation (see [2])
%     shiftAxis             shift along this axis ('x','y' (default),'z')
%     shiftOffset           shift ears (L - R) by these values 
%                           (default:  [-0.0725 0.0725])
%
% OUTPUT:
%     itaHRTF object
%     .freqData: interpolated / range-extrapolated HRTFs for defined field points
%     .timeData: interpolated / range-extrapolated HRIRs for defined field points
%     .dirCoord: itaCoordinates object
%
% 
%
% [1] Pollow, Martin et al., "Calculation of Head-Related Transfer Functions
%     for Arbitrary Field Points Using Spherical Harmonics Decomposition",
%     Acta Acustica united with Acustica, Volume 98, Number 1, January/February 2012,
%     pp. 72-82(11)
%
% [2] Richter, Jan-Gerrit et al. "Spherical harmonics based hrtf datasets: 
%     Implementation and evaluation for real-time auralization",
%     Acta Acustica united with Acustica, Volume 100, Number 4, July/August 2014,
%     pp. 667-675(9)
%
% Author:  Florian Pausch <fpa@akustik.rwth-aachen.de>
% Version: 2016-02-05

sArgs           = struct('order',50,'epsilon',1e-8,'shiftToEar',false,'shiftAxis','y','shiftOffset',[-0.0725 0.0725]);
sArgs           = ita_parse_arguments(sArgs,varargin,2);
if isempty(varargin) || ~isa(varargin{1},'itaCoordinates')
    error('itaHRTF:interp', ' An itaCoordinate object is needed!')
end
field_in        = varargin{1};

% only take unique direction coordinates (round to 0.01deg resolution) 
tempfield       = unique(round([field_in.phi_deg*100 field_in.theta_deg*100]),'rows'); % may cause problems with older Matlab versions (<=R2013)!
tempfield = tempfield./100;
temp_r          = field_in.r(1);
field           = itaCoordinates(size(tempfield,1));
field.r         = repmat(temp_r,size(tempfield,1),1);
field.phi_deg   = tempfield(:,1);
field.theta_deg = tempfield(:,2);

N               = sArgs.order;
epsilon         = sArgs.epsilon;                                   % regularization parameter
k               = this.wavenumber;                             % wave number
k(1)            = eps;
% add eps to avoid NaN's
% Nmax            = floor(sqrt(this.nDirections/4)-1);
Nmax            = ceil(sqrt(this.nDirections/2)-1);
% Nmax            = 50;

if N>Nmax
   N=Nmax;
   warning(['Order of SH reconstruction matrix was set to N = Nmax = ',num2str(Nmax),'.'])
end

% construct vector of length (N+1)^2 regularization weights and,
% if needed, spherical hankel functions of second kind (for r0 and r1)
if ~isequal(this.dirCoord.r(1),field.r(1))
    kr0 = k*this.dirCoord.r(1);                % measurement radius [m]
    kr1 = k*field.r(1);                        % extrapolation radius [m]
    
    hankel_r0 = ita_sph_besselh(ita_sph_linear2degreeorder(1:Nmax),2,kr0);
    hankel_r1 = ita_sph_besselh(ita_sph_linear2degreeorder(1:Nmax),2,kr1);
    hankel_div = hankel_r1 ./ hankel_r0;
    
    hankel_rep = hankel_div(:,1);
end

nSH = (Nmax+1).^2;
I = sparse(eye(nSH));
n = ita_sph_linear2degreeorder(1:round(nSH)).';

%% move data from earcenter
copyData = this;
if sArgs.shiftToEar
    ear_d       =   sArgs.shiftOffset;
    for ear=1:2
        movedData = moveFullDataSet(this.ch(ear:2:this.nChannels),sArgs,ear_d(ear),1);
        copyData.freqData(:,ear:2:copyData.nChannels) = movedData;
    end
end

%% Weights
[~,w]= this.dirCoord.spherical_voronoi;         % calculate weighting coefficients (Voronoi surfaces <-> measurement points)

W = sparse(diag(w));                                      % diagonal matrix containing weights
D = I .* diag(1 + n.*(n+1));                              % decomposition order-dependent Tikhonov regularization
Y = ita_sph_base(this.dirCoord,Nmax,'orthonormal',true);  % calculate real-valued SHs using the measurement grid

%% Calculate HRTF data for field points
if Nmax > 25
    ita_disp('[itaHRTF.interp] Be patient...')
end
    
% init.
hrtf_arbi = zeros(this.nBins,2*field.nPoints); % columns: LRLRLR...
for ear=1:2
    
    % SH transformation
    freqData_temp   = copyData.freqData(:,ear:2:end);
    a0              = (Y.'*W*Y + epsilon*D) \ Y.'*W * freqData_temp.';

%     %% test the sh transformation results by plotting both spatial and sh data with surf
%     s = itaSamplingSph(field);
%     s.nmax = Nmax;
%     figure
%     surf(s,a0(:,10))
%     figure
%     surf(s,freqData_temp(10,:))
    
    % range extrapolation
    if ~isequal(this.dirCoord.r(1),field.r(1))
        % calculate range-extrapolated HRTFs
        a1 = a0 .* hankel_rep.';
        
        %%% test here to see extrapolation results in spatial domain
        %     surf(s,a1(:,10))
        
        % reconstruction to spatial data
        Yest = ita_sph_base(field,N,'orthonormal',true);  % use real-valued SH's
        hrtf_arbi(:,ear:2:end) = (Yest*a1(1:(N+1)^2,:)).';             % interpolated + range-extrapolated HRTFs
    else
        % reconstruction to spatial data
        Yest = ita_sph_base(field,N,'orthonormal',true);  % use real-valued SH's
        hrtf_arbi(:,ear:2:end) = (Yest*a0(1:(N+1)^2,:)).';             % interpolated HRTFs
    end
end


% set new direction coordinates
sph                         = zeros(field.nPoints*2 ,3);
sph(1:2:end,:)              = field.sph;
sph(2:2:end,:)              = field.sph;

% write new HRTF data set
cThis                       = this;
cThis.freqData = hrtf_arbi;
cThis.channelCoordinates.sph= sph;

% channelnames coordinates
cThis.channelNames = ita_sprintf('%s (\\theta=%2.0f, \\phi=%2.0f)',...
    repmat(['L'; 'R'],cThis.dirCoord.nPoints, 1),...
    cThis.channelCoordinates.theta_deg,...
    cThis.channelCoordinates.phi_deg);

%% move back to head center
if sArgs.shiftToEar
    ear_d_back       =  ear_d;
%     movedData = zeros(size(hrtf_arbi));
    for ear=1:2
        movedData = moveFullDataSet(cThis.ch(ear:2:cThis.nChannels),sArgs,ear_d_back(ear),2);
        cThis.freqData(:,ear:2:cThis.nChannels) = movedData;
    end    
end

if ~isequal(cThis.dirCoord.r(1),field.r(1))%???
    cThis.dirCoord.r = field.r;
end

if N > 25
    ita_disp('[itaHRTF.interp] ...calculation finished!')
end
end



function [ data ] = moveFullDataSet(data,options,offsetShift,mode)
    fullCoords = data.channelCoordinates;
    freqVector = data.freqVector;
    shiftedData = zeros(size(data.freqData));
    axis = options.shiftAxis;
    for index = 1:length(freqVector)
        shiftedData(index,:) = moveHRTF(fullCoords,data.freqData(index,:),freqVector(index),axis,offsetShift,mode);
    end


    data = shiftedData;
end


function [data] = moveHRTF(s, data, frequency, axis, offset,mode)
    % the offset is given in m
    
    origAxis = s.r;
    
    if (size(data,2) > size(data,1))
       data = data.'; 
    end
    offset = real(offset); % ??
    switch axis
        case 'x' 
            s.x = s.x + offset;
        case 'y'
            s.y = s.y + offset;
        case 'z'
            s.z = s.z + offset;
    end
    
    newAxis = s.r;
    k = 2*pi*frequency/340;
    % the phase is moved by the difference of the axis points
    switch mode
        case 1
            data = data .* exp(1i*k*(newAxis - origAxis));
        case 2
            data = data .* exp(1i*k*(origAxis - newAxis));     
    end
    % amplitude manipulation did not yield better results
%     data = data .* newAxis ./ origAxis;

end

function interpolatedHRTF = miinterpolateHRTF(hrtfData,sourcePosition,...
    desiredPosition,varargin)

%#codegen

    % Check if hrtfData, sourcePosition, and
    % desiredPosition are provided
    narginchk(3,5);
    
    % Parse the input after validating
    [hrtfInterpolationObject,interpolationAlgorithm] = parseInput(hrtfData,...
        sourcePosition,desiredPosition,varargin);
    
    % Pre-initialize the variables
    numOfDesiredPos = size(desiredPosition);
    sizeOfData = size(hrtfData);
    interpolatedHRTF = zeros(numOfDesiredPos(1), sizeOfData(2),...
        sizeOfData(3),class(hrtfData));
    
    % Pass a pair of Az, El from vector inputs
    for idx = 1:numOfDesiredPos(1)
        % Codegen requires the variable's size to be defined before assigning.
        % Hence to avoid for-loops using an extra variable
        interimResult = performInterpolation(hrtfInterpolationObject,...
            hrtfData, desiredPosition(idx,:), interpolationAlgorithm);
        
        interpolatedHRTF(idx,1:sizeOfData(2),1:sizeOfData(3)) = ...
            interimResult(1:sizeOfData(2),1:sizeOfData(3));
    end
    
end

function [hrtfInterpolationObject, algorithm] = parseInput(hrtfData,sourcePosition,desiredPosition,argIn)
    
    % Set the puDomain property
    % If HRTFData is real, domain is time. Else frequency
    if(isreal(hrtfData))
        domain = 'time';
    else
        domain = 'frequency';
    end

    % Default algorithm is bilinear
    algorithm = 'bilinear';
    
    if(~isempty(argIn))       
        % Check if only one name-value pair is provided
        validateattributes(argIn,{'cell'},{'numel',2},...
            'interpolateHRTF','Name-Value pair')
        
        % Check if the Name in Name-Value pair is appropriate
        validatestring(argIn{1},{'Algorithm'},'interpolateHRTF','Name');
               
        % Validate the contents of 'Value'
        algorithm = validatestring(argIn{2},...
            {'vbap','bilinear'},'interpolateHRTF',...
            'the value of ''Algorithm''');        
    end
    
%     hrtfInterpolationObject = HRTFInterpolation(...
%         'InterpolationAlgorithm',algorithm ,...
%         'SourcePosition',sourcePosition,...
%         'puDomain',domain);
    hrtfInterpolationObject.InterpolationAlgorithm = algorithm;
    hrtfInterpolationObject.SourcePosition = setSourcePos(sourcePosition);
    hrtfInterpolationObject.pDomain = domain;
end


function SourcePosition = setSourcePos(val)
        % Internally, always use 0 to 360 convention for source positions
        val(:,1) = mod(val(:,1), 360);
        SourcePosition = val;
end        

    
%%
        
        function interpolatedHRTF = performInterpolation(obj,hrtfData,...
                                  desiredPosition,interpolationAlgorithm)
            % Set the properties of the interpolateHRTF class
            % Internally, always use 0 to 360 convention for desired position
            desiredPosition(1,1) = mod(desiredPosition(1,1), 360);
            desiredAz = desiredPosition(1,1);
            desiredEl = desiredPosition(1,2);                                     
            
            % Call the appropriate interpolation function
            if(strcmpi(interpolationAlgorithm,'bilinear'))
                % Find the 3 nearest reference sources for the desired
                % Az and El for Bilinear
                obj = getReferenceSourceBL(obj,desiredPosition);
                interpolatedHRTF = performBLInterpolation(obj,hrtfData,desiredAz,desiredEl);
            
            elseif(strcmpi(interpolationAlgorithm,'vbap'))
                % Find the 3 nearest reference sources for the desired
                % Az and El for VBAP
                obj = getReferenceSourceVBAP(obj,desiredPosition);
                interpolatedHRTF = performVBAPInterpolation(obj,hrtfData,desiredAz,desiredEl);
            else
                interpolatedHRTF = ...
                    zeros(size(hrtfData), class(hrtfData));
            end            
        end        
    
    
        function outVal = getCartesian(azimuth, elevation,r)
            % Az, el are in degrees. So convert to radians to use the sph2cart function
            [x ,y ,z] = sph2cart(azimuth * pi/180, elevation * pi/180,r);
            outVal = [x y z];
        end
        
        function callError (interpolationAlgorithm)
         coder.internal.error('audio:interpolateHRTF:NoReferenceSources',...
               interpolationAlgorithm); 
        end
        
    
    
%% Interpolação BILINEAR       
function obj = getReferenceSourceBL(obj,desiredPosition)            
    % Get the nearest loud speaker Az and El from the desired Az
    % and El            

    % Get the difference between source pos and the desired pos
    N = length(obj.SourcePosition);
    diffSource = zeros(N,2,class(obj.SourcePosition));

    diffSource(1:N,1) = (obj.SourcePosition(1:N,1) - desiredPosition(1,1));
    diffSource(1:N,2) = (obj.SourcePosition(1:N,2) - desiredPosition(1,2));            

    [~, sortedID] = unique(abs(diffSource),'rows');

    % Point A
    obj.pRefSources.Az(1) = obj.SourcePosition(sortedID(1),1);
    obj.pRefSources.El(1) = obj.SourcePosition(sortedID(1),2);
    obj.pRefSources.Index(1) = sortedID(1);
    count = 2;

    sAz = obj.SourcePosition(sortedID,1);
    sEl = obj.SourcePosition(sortedID,2);

    % Point B            
    % Point A and B should have same El (are you sure?)
    id = ( (obj.pRefSources.Az(1) ~= sAz));
    idB = find(id);

    if(isempty(idB))
        % If the appropriate references sources cannot be found,
        % throw an error
        obj.callError('Bilinear');
    end

    obj.pRefSources.Az(2) = obj.SourcePosition(sortedID(idB(1)),1);
    obj.pRefSources.El(2) = obj.SourcePosition(sortedID(idB(1)),2);
    obj.pRefSources.Index(2) = sortedID(idB(1));

    % Point C
    id = find(obj.pRefSources.El(1:count-1) ~= sEl);
    zz = find(id > idB(1));
    idC = id(zz(1));

    if(isempty(idC))
        % If the appropriate references sources cannot be found,
        % throw an error
        obj.callError('Bilinear');
    end

    obj.pRefSources.Az(3) = obj.SourcePosition(sortedID(idC(1)),1);
    obj.pRefSources.El(3) = obj.SourcePosition(sortedID(idC(1)),2);
    obj.pRefSources.Index(3) = sortedID(idC(1));    
end
        
        
function interpolatedIR = performBLInterpolation(obj, hrtfData, desiredAz, desiredEl)
    % Get the 3-point bilinear interpolation weights
    % phi - Elevation
    % theta - Azimuth

    % Point A and B are assumed to have same Elevation

    % wC = [phi(p) - phi(a) / phi(c) - phi(a)]
    weightC = (desiredEl - obj.pRefSources.El(1)) / (obj.pRefSources.El(3) - obj.pRefSources.El(1));

    % wB = [(theta(p) - theta(a)) - (wC*(theta(c) - theta(a)))] / [theta(b) - theta(a)]
    weightB = ( (desiredAz - obj.pRefSources.Az(1)) - ( weightC * (obj.pRefSources.Az(3)- obj.pRefSources.Az(1) ) ))/ (obj.pRefSources.Az(2) - obj.pRefSources.Az(1));

    % wA = 1 - wB - wC
    weightA = 1 - weightB - weightC;

    refSourceLength = 3;
    nSamples = size(hrtfData,3);

    % Pre-initializing values to zeros
    hLeftMin = zeros(refSourceLength,nSamples,class(hrtfData));
    hRightMin = zeros(refSourceLength,nSamples,class(hrtfData));
    interpolatedIR = zeros(2,nSamples,class(hrtfData));

    % Output is same type as the input
    hLeftDelay = zeros(refSourceLength,1,class(hrtfData));
    hRightDelay = zeros(refSourceLength,1,class(hrtfData));

    for ix = 1: refSourceLength                
        % Find the index of the known HRTF 
        idx = find(obj.SourcePosition(:,1,1) == obj.pRefSources.Az(ix)...
            & obj.SourcePosition(:,2,1) == obj.pRefSources.El(ix));
        hLeft = squeeze(hrtfData(idx, 1, :)); % Get the Left ear HRIR
        hRight = squeeze(hrtfData(idx, 2, :)); % Get the Right ear HRIR

        % Get the Minimum Phase of the HRIR
        % In order to avoid undesired destructive interference, the interpolation
        % procedure given by equation above should be performed on
        % the minimum-phase versions of the measured HRIRs

        [hL , hLeftDelay(ix)]= getMinimumPhase(obj,hLeft);
        [hR, hRightDelay(ix)] = getMinimumPhase(obj,hRight);

        % Codegen requires the variable size before assigning.
        % Hence to avoid for-loops using an extra variable
        hLeftMin(ix,1:nSamples) = hL(1:nSamples,1);
        hRightMin(ix,1:nSamples)= hR (1:nSamples,1);

    end

    % Perform the interpolation on the Minimum phase HRIR
    LinterpolatedIR = weightA*hLeftMin(1,:) + weightB*hLeftMin(2,:) + weightC*hLeftMin(3,:);
    RinterpolatedIR = weightA*hRightMin(1,:) + weightB*hRightMin(2,:) + weightC*hRightMin(3,:);

    % Get the interpolated delay
    leftDelay = floor( weightA*hLeftDelay(1,:) + weightB*hLeftDelay(2,:) + weightC*hLeftDelay(3,:));
    rightDelay = floor(weightA*hRightDelay(1,:) + weightB*hRightDelay(2,:) + weightC*hRightDelay(3,:));

    if(strcmpi(obj.pDomain , 'time'))
        % Time-Domain - compute the actual Interpolated HRIR using delay
        % This will be convolution of hmin[k] with a unit impulse
        % delayed by tau in samples

        if(leftDelay >= 0)
            interpolatedIR(1,leftDelay+1:end) = ...
                LinterpolatedIR(1,1:end-leftDelay);
        else
            lDelay = abs(leftDelay);
            interpolatedIR(1,1:end-lDelay+1) = ...
                LinterpolatedIR(1,lDelay:end);
        end

        if (rightDelay >= 0)
            interpolatedIR(2,rightDelay+1:end) = ...
                RinterpolatedIR(1,1:end-rightDelay);
        else
            rDelay = abs(rightDelay);
            interpolatedIR(2,1:end-rDelay+1) = ...
                RinterpolatedIR(1,rDelay:end);
        end                                

    elseif(strcmpi(obj.pDomain , 'frequency'))
        % Frequency Domain - - compute the actual Interpolated HRTF using delay
        % Get the delay constant: exp(-j*2*pi*k*Delay/N)
        % Multiply this constant with LinterpolatedIR and
        % RinterpolatedIR
        N = length(LinterpolatedIR);
        k = 0:N-1;
        LdelayConstant = exp(-1i*2*pi*k*leftDelay / N);
        RdelayConstant = exp(-1i*2*pi*k*rightDelay / N);

        interpolatedIR(1,:) = LdelayConstant .* LinterpolatedIR;
        interpolatedIR(2,:) = RdelayConstant .* RinterpolatedIR;
    end
end                               

function [hMin, tau] = getMinimumPhase(obj,hInput)
    % To find the minimum phase HRIR using Hilbert transform
    % The minimum phase response and magnitude spectrum are connected
    % via the Hilbert transform
    % phaseRespMin = imaginary(HilbertTransform(-ln(magnitude of hSignal)))

    hSignal = hInput;

    if(strcmpi(obj.pDomain , 'time'))
        % If time-domain, we need to compute the FFT of the HRIR
        hSignal = fft(hInput);
    end

    hMag = abs(hSignal);            
    hLnMag = log(hMag);            
    hHilbert = hilbert(-hLnMag);            
    hMinPhase = imag(hHilbert);

    hMin = ifft((hMag.*exp(1i*hMinPhase)));
    hMin = real(hMin);


    % To find the delay tau
    % For the estimation of the absolute delay tau , using the spectra
    % of an impulse response (H[]) and its minimum phase version Hmin[] for
    % calculating the cross-correlation
    hMinFFT = fft(hMin);

    [~, tau] = max(ifft(hSignal.*conj(hMinFFT)));

    if(strcmpi(obj.pDomain , 'frequency'))
        % For Frequency Domain, we need the FFT of the minimum
        % phase filter
        hMin = hMinFFT;
    end            
end
        

%% Interpolação VBAP
function obj = getReferenceSourceVBAP(obj, desiredPosition)
    % Get the nearest loud speaker Az, and El from the desired Az
    % and El            

    N = length(obj.SourcePosition);
    diffSource = zeros(N,2,class(obj.SourcePosition));

    diffSource(1:N,1) = (obj.SourcePosition(1:N,1) - desiredPosition(1,1));
    diffSource(1:N,2) = (obj.SourcePosition(1:N,2) - desiredPosition(1,2));

    [~, sortedID] = unique(abs(diffSource),'rows');

    % Point A
    obj.pRefSources.Az(1) = obj.SourcePosition(sortedID(1),1);
    obj.pRefSources.El(1) = obj.SourcePosition(sortedID(1),2);
    obj.pRefSources.Index(1) = sortedID(1);

    sAz = obj.SourcePosition(sortedID,1);
    sEl = obj.SourcePosition(sortedID,2);

    % Point B            
    zB = (obj.pRefSources.Az(1) == sAz) |...
        (obj.pRefSources.El(1) == sEl);

    idx = find(zB == 0);

    % This case is unlikely, but if the 3 nearest reference
    % sources are not found for VBAP, we cannot proceed ahead
    if(isempty(idx))
        obj.callError('VBAP');
    end

    obj.pRefSources.Az(2) = obj.SourcePosition(sortedID(idx(1)),1);
    obj.pRefSources.El(2) = obj.SourcePosition(sortedID(idx(1)),2);
    obj.pRefSources.Index(2) = sortedID(idx(1));

    % Point C
    zC = (obj.pRefSources.Az(2) == sAz) | (obj.pRefSources.El(2) == sEl);
    idX = find(zB == 0 & zC == 0);

    if(isempty(idX))
        obj.callError('VBAP');
    end

    obj.pRefSources.Az(3) = obj.SourcePosition(sortedID(idX(1)),1);
    obj.pRefSources.El(3) = obj.SourcePosition(sortedID(idX(1)),2);
    obj.pRefSources.Index(3) = sortedID(idX(1));                                              

end

function interpolatedIR =  performVBAPInterpolation(obj,hrtfData,desiredAz,desiredEl)
    % Pre-initializing values to zeros
    refSourceLength = 3;
    nSamples = size(hrtfData,3);

    hLeft = zeros(refSourceLength,nSamples,class(hrtfData));
    hRight = zeros(refSourceLength,nSamples,class(hrtfData));

    interpolatedIR = zeros(2,nSamples,class(hrtfData));


    % Value of radius does not matter since the radius does not change
    % and same value is used to compute p and L matrices
    radius = 1;

    % Compute p
    p = getCartesian(desiredAz, desiredEl,radius);

    % Compute l

    l1 = getCartesian(obj.pRefSources.Az(1), obj.pRefSources.El(1),radius);
    l2 = getCartesian(obj.pRefSources.Az(2), obj.pRefSources.El(2),radius);
    l3 = getCartesian(obj.pRefSources.Az(3), obj.pRefSources.El(3),radius);

    l = [l1;l2;l3];

    % Calculate the VBAP Gain
    % g = p^T * l^-1
    % gL = p -> g = mrdivide(p,L) or g = p/l;                                        
%     g = p/l;     % FUNÇÂO ORIGINAL 
    g = p*pinv(l); % ALTERAÇÂO EM RELAÇÂO A FUNÇÂO ORIGINAL 
    
    % Normalize the gains
    g = g/(sqrt(sum(g.^2)));            

    % mrdivide can result in unreliable results when A is 
    % singular. We could use lsqminnorm instead but it is not
    % supported for codegen. Hence, suppressing the 'MATLAB:nearlySingularMatrix'
    % warning since we are validating if the results are acceptable
    % later.

    coder.extrinsic('lastwarn', 'warning');
    [~,msgId] = lastwarn;
    warnStruct = warning('off',msgId);

    % Validate that gain 'g' is a finite value and not-NAN
    if(isnan(g) | isinf(g)) %#ok<OR2>           
        % While selecting the sources, all care has been taken to choose
        % sources which does not result in singular matrices. Even then 
        % the sources result in singular matrices, we need different data
        obj.callError('VBAP');
    end

    % Restore previous warning state
    warning(warnStruct);            

    % Get the known HRIR's
    for ix = 1:refSourceLength
        % Find the index of the known HRTF 
        idx = find(obj.SourcePosition(:,1,1) == obj.pRefSources.Az(ix)...
            & obj.SourcePosition(:,2,1) == obj.pRefSources.El(ix));

        % Typically, idx will be a scalar. But, if the dataset has
        % repeated measurements for the same Az and El, idx will be
        % a matrix. Hence, accept only the first value in idx
        hL = squeeze(hrtfData(idx(1), 1, :)); % Get the Left ear HRIR;
        hR = squeeze(hrtfData(idx(1), 2, :)); % Get the Right ear HRIR

        hLeft(ix,1:nSamples) = hL(1:nSamples,1);
        hRight(ix,1:nSamples) = hR(1:nSamples,1);

    end            
    interpolatedIR(1,:) = g * hLeft;
    interpolatedIR(2,:) = g * hRight;            
end


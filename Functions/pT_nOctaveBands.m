function [bandCenter,bandValue] = pT_nOctaveBands(inputFreq,ValueData,varargin)
%pT_nOctaveBand Provides midband frequency values and calculated data
%values on nth octave bands distributed in compliance with ANSI.S1.11.2004.
%Variable distances in the frequency vector are supported. inputFreq vector
%must be ascending in order. No conversions (for instance, from/to Pa 
%to/from dB) are made by the code in order to provide enhanced control 
%for the user. Please notice that, unlike traditional narrow-to-band 
%conversion, this algorithm splits the values between band center
%frequencies among both bands, applying a window of the form 
%cos(0.5*n*pi)^2, where "n" is the base 2 logarithmic displacement 
%from the band center frequency. The trigonometric property encompassed 
%ensures that no losses are possible through this method. The input value is then 
%multiplied by the given window coefficient to equilibrate its contribution
%between frequencies. On averaging mode, this is done by the integral of 
%the window function on the input value corresponding bandwidth; On summing mode,
%this is done by the sheer multiplication of the window function value on 
%the frequency value given.
%
%       INPUT PARAMETERS
% inputFreq -> Nx1 matrix of center frequencies from input data
% ValueData -> Nx1 matrix of values (Pa, W, W/m^2...) corresponding the
% previously specified frequencies
%
% By default the output will be given in
% b: 3 - Third Octave Bands
% Mode: 'avg' - Averages the values within each nTh octave
% band
% plots: false - No automatic plots will be generated.
% 
%       ADDITIONAL PARAMETERS
% b (integer and positive) ->  Bandwidth designator, reciprocal of a 
% positive integer to designate the fraction of an octave band. A value of
% 3 provides the result in third octave bands, a value of 5 provides the 
% result in fifth octave bands, and so on...
%
% Mode ('avg' or 'sum') -> Defines the mode of evaluation. 'avg' averages 
% the values within the band, while 'sum' considers that each frequency 
% component contributes individually for the band's overall value. Average 
% mode may be useful if, for instance, post-processing over simulation data 
% of a component's transmission loss is being done. Sum mode may be more 
% useful to applications such as measurements, where narrow band data is
% introduced and its values in nth octave bands are, indeed, meant to 
% contribute with one another.
%
% plots (false or true) -> Plots processing results in a separate
% window if logical true is specified.
%
% ex: pT_nOctaveBand([200 600 2000], [3.1546 2.5165 9.854]) --Puts the
% first vector (frequencies) and second vector (Pa/W/Wm^-2 values) in third
% octave bands form.
% ex: pT_nOctaveBand((20:2:20000),rand(1,9991),'b',10) -- Ensures that the
% output will be given 10th octave bands.
% ex: pT_nOctaveBand((20:2:20000),rand(1,9991),'b',1,'plots',true) -- Provides
% the same type of outputs as above in whole octave bands and plots the 
% results.
% ex: pT_nOctaveBand((20:2:20000),rand(1,9991),'Mode','sum') - Provides the
% result summing each individual component contribution for each band,
% instead of averaging their values amongst them.
%
% FULL SYNTAX: [FreqResults,ValueResults] =
% pT_nOctaveBand(FreqArray,ValueArray,'b',double,'Mode',string,'plots',boolean);
% where:
% double - Any integer value greater than or equal to 1
% string - 'avg' or 'sum' allowed only
% boolean - true or false values, without quotes.
%
% IMPORTANT: For realtime tracking, do not exceed more than 500 frequency
% points per sample, otherwise processing times may outrun your sampling
% rate. An average of 0.1sec per sample is expected when working with 500
% data points. Considering using a separate vector, with fewer data points,
% if your input may exceed those values.

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%INPUT VALUES FOR TESTING - START
% b = 3;
% inputFreq = [20 40 60 80 100 120 140 160 180 200 220 240 260 280 300 320 340 360 380 400 420 440 460 480 500 520 540 560 580 600 620 640 660 680 700 720 740 760 780 800 820 840 860 880 900 920 940 960 980 1000 1020 1040 1060 1080 1100 1120 1140 1160 1180 1200 1220 1240 1260 1280 1300 1320 1340 1360 1380 1400 1420 1440 1460 1480 1500 1520 1540 1560 1580 1600 1620 1640 1660 1680 1700 1720 1740 1760 1780 1800 1820 1840 1860 1880 1900 1920 1940 1960 1980 2000 2020 2040 2060 2080 2100 2120 2140 2160 2180 2200 2220 2240 2260 2280 2300 2320 2340 2360 2380 2400 2420 2440 2460 2480 2500 2520 2540 2560 2580 2600 2620 2640 2660 2680 2700 2720 2740 2760 2780 2800 2820 2840 2860 2880 2900 2920 2940 2960 2980 3000];
% ValueData = [2.15787154351975 5.50910336776418 8.24810621910948 10.3852765393456 12.0880928069250 13.4764727879142 14.6282351492306 15.5947324466101 16.4109897059518 17.1016592301036 17.6845490302261 18.1727771168130 18.5761296147740 18.9019414786028 19.1556801433281 19.3413370829109 19.4616898757354 19.5184726049571 19.5124772236968 19.4435985450564 19.3108282716027 19.1121973486084 18.8446596482629 18.5039022612325 18.0840567206945 17.5772684711860 16.9730537920353 16.2573247943115 15.4108758476928 14.4069631620977 13.2073054096862 11.7552842479692 9.96439720571211 7.70154945731211 4.79474745181861 1.40342467784380 0.153550040450445 2.99376719327999 6.25874958770432 8.84956412327973 10.8766891147885 12.5010676805012 13.8318273387441 14.9396416720048 15.8714412792712 16.6594630886054 17.3265384861789 17.8892411833143 18.3598219293658 18.7485599580903 19.0600821724185 19.3005831431679 19.4737533999101 19.5821156280716 19.6271774980075 19.6095190638370 19.5288255386930 19.3838693699504 19.1724394578529 18.8912088367556 18.5355237835582 18.0990851579850 17.5734735369987 16.9474374678725 16.2058077126217 15.3277976868370 14.2842576681772 13.0330841328486 11.5113238257373 9.62172037005920 7.21451509029879 4.11612489126520 0.763166396820580 0.599392504807653 3.89982411081434 7.05844390687688 9.51965293275296 11.4510621163089 13.0063584674922 14.2857746412979 15.3540839552287 16.2545031774786 17.0168870895033 17.6624870411818 18.2067907151142 18.6612758258433 19.0345282203503 19.3329736839825 19.5613657165149 19.7231128256268 19.8204954418503 19.8548026580181 19.8264064694558 19.7347826322068 19.5784806241358 19.3550390674372 19.0608360870328 18.6908549328648 18.2383315346856 17.6942286382679 17.0464437037684 16.2785911254610 15.3680761249162 14.2829420291693 12.9765146420880 11.3780210309885 9.37636695115089 6.79926694881442 3.47424455560014 0.280468835760290 1.31408273063171 4.89660585915041 7.94580886592248 10.2946259049822 12.1435845708140 13.6389766554364 14.8736847270383 15.9076017158816 16.7800967829245 17.5209845666808 18.1482707502326 18.6770653132536 19.1181693340052 19.4796605737720 19.7675727882702 19.9863460589634 20.1391280163884 20.2279718261865 20.2539582594459 20.2172575136454 20.1171383185558 19.9519252463672 19.7188981639060 19.4141121056984 19.0325541979820 18.5660171506031 18.0055928761450 17.3380247159350 16.5455482880198 15.6034869578951 14.4765715909460 13.1125564552175 11.4307066608382 9.30113113612342 6.51755705936837 2.90268143370247 0.0228967718824938 2.29433707854459 6.04763070374218 9.00185684246989];
% warning("Test data in use for nth octave band function. Input arguments are not being used.")
%INPUT VALUES FOR TESTING - END
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

for i = 1:length(inputFreq)
    FreqData(i).Value = inputFreq(i);
end

b = 3; %Third octave bands - Default value
Mode = 'avg'; %Averagin mode - Default Value
plots = false;

%Checking additional arguments available
if length(varargin) > 0
    for i = 1:length(varargin)
        if length(varargin{i}) == length('b')
            if varargin{i} == 'b'
                b = varargin{i+1};
            end
        end
        if length(varargin{i}) == length('Mode')
            if eq(varargin{i},'Mode')
                Mode = varargin{i+1};
            end
        end
        if length(varargin{i}) == length('plots')
            if varargin{i} == 'plots'
                plots = varargin{i+1};
            end
        end
    end
end

%Checking for an integer input of 'b'
if ne(rem(b,1),0)
    error("Only integer fractions of octave bands are allowed. Please recheck your inputs. Octave band fraction specified: 1/%f", b);
end

%Checking for correct string at Mode
if ne(Mode,'avg') & ne(Mode,'sum')
    error("Only strings 'sum' and 'avg' are allowed in Mode parameter");
end

%Checking for correct values for plots
if ne(plots,true) & ne(plots,false)
    error("Only booleans true or false are accepted at plots parameter");
end


G = 10^(3/10);   % ANSI.S1.11.2004 - Octave ratio G (base 10) frequency ratio
f_ref = 1000;       % ANSI.S1.11.2004 - Standard reference frequency [Hz]

%% Main calculations
%Calculating midband frequencies for odd b's
if rem(b,2)
    for i = 1:2000
        if (G^(-1/(2*b))*f_ref*(G^((i-30)/b)) > 20000) %That is: if the lower edge of the frequency band to be calculated is higher than 20000Hz, stop calculating more bands
            break
        end
        band(i).center = f_ref*(G^((i-30)/b));
    end
end

%Calculating midband frequencies for even b's
if (rem(b,2)+1)
    for i = 1:2000
        if (G^(-1/(2*b))*f_ref*(G^((i-30)/b)) > 20000) %That is: if the lower edge of the frequency band to be calculated is higher than 20000Hz, stop calculating more bands
            break
        end
        band(i).center = f_ref*(G^((2*i-59)/(2*b)));
    end
end

%% Imposing nominal midband frequencies according to table A1 of the 
%  standard when b is equal to 1 or 3.
if b == 3
    clear band
    nBf = [25 31.5 40 50 63 80 100 125 160 200 250 315 400 500 630 800 1000 1250 1600 2000 2500 3150 4000 5000 6300 8000 10000 12500 16000 20000];
    for i = 1:length(nBf)
        band(i).center = nBf(i);
    end
    
end

if b == 1
    clear band
    nBf = [31.5 63 125 250 500 1000 2000 4000 8000 16000];
    for i = 1:length(nBf)
        band(i).center = nBf(i);
    end
end
    
numberofbands = length(band);
    
%Calculating edge frequencies
for i = 1:numberofbands
    band(i).edge1 = G^(-1/(2*b))*(band(i).center);  %Lower limit
    band(i).edge2 = G^(1/(2*b))*(band(i).center);   %Higher limit
end

%%  Nominal midband frequencies for bandwidth designator 1/4 to 1/24 inclusive. Higher precision for b>24 also included.
%   This section rounds midband freq values according to the standard
%   The standard only mentions that for b > 24, a higher precision is
%   demanded, but this higher precision is not specified. Therefore,
%   additional 3 digits will be included.
if ge(b,4) 
    for i = 1:numberofbands
        %Getting digit to the left
        leftmost = (num2str(band(i).center));
        leftmost = str2num(leftmost(1));

        %Rounding rule for leftmost digit between 1-4 (inclusive)
        if ge(leftmost,1) && le(leftmost,4)
            digits = num2str(band(i).center);
            decDigitPos = length(digits)+1;    %Assuming there's no decimal digit in the string, this provides that its initial position to be considered is one character beyond the num2str result
            for j = 1:length(digits)            %Finding and tagging decimal point position
                if digits(j) == '.'             %Finding and tagging decimal point position
                    decDigitPos = j;            %Finding and tagging decimal point position
                end
            end

            %Rounding midband frequencies
            %Compliance with ANSI.S1.11.2004 for rounding midband frequencies for 
            %bandwidth designators 1/4 to 1/24 inclusive
            n = 4 - decDigitPos; %number of decimal digits to round up to (or integers to lose significance)
           
            if b > 24
                n = n + 3; % Increasing precision for b>24
            end
            
            band(i).center = round(band(i).center,n);
        end

        %Rounding rule for leftmost digit between 5-9 (inclusive)
        if ge(leftmost,5) && le(leftmost,9)
            digits = num2str(band(i).center);
            decDigitPos = length(digits)+1;    %Assuming there's no decimal digit in the string, this provides that its initial position to be considered is one character beyond the num2str result
            for j = 1:length(digits)            %Finding and tagging decimal point position
                if digits(j) == '.'             %Finding and tagging decimal point position
                    decDigitPos = j;            %Finding and tagging decimal point position
                end
            end

            %Rounding midband frequencies
            %Compliance with ANSI.S1.11.2004 for rounding midband frequencies for 
            %bandwidth designators 1/4 to 1/24 inclusive
            n = 3 - decDigitPos; %number of decimal digits to round up to (or integers to lose significance)
            
            if b > 24
                n = n + 3; % Increasing precision for b>24
            end
            
            band(i).center = round(band(i).center,n);
        end
    end
end

%% Evaluation edges for first and last bands - A different band edge for mathematical purposes
evaluationedge1 = zeros(1,1);
evaluationedge2 = zeros(1,1);
evaluationedge1 = evaluationedge1 + (band(1).center * 2^(-log2((band(2).center)/band(1).center))); 
evaluationedge2 = evaluationedge2 + (band(length(band)).center * 2^(-log2((band(length(band)-1).center)/band(length(band)).center))); 


%% Main calculations (2) - Establishing the start and end of the "rectangles" occupied by input frequency data

for i = 1:length(FreqData)
    if i > 1 && i < length(FreqData)
        FreqData(i).edge1 = FreqData(i).Value * 2^(0.5*log2((FreqData(i-1).Value)/FreqData(i).Value));
        FreqData(i).edge2 = FreqData(i).Value * 2^(0.5*log2((FreqData(i+1).Value)/FreqData(i).Value));
    end
    
    if i==1
        FreqData(i).edge1 = FreqData(i).Value * 2^(-0.5*log2((FreqData(i+1).Value)/FreqData(i).Value)); %This assures the lower half of the first rectangle has the same size of the upper half (which has data to determine it)
        FreqData(i).edge2 = FreqData(i).Value * 2^(0.5*log2((FreqData(i+1).Value)/FreqData(i).Value));
    end
    
    if i==length(FreqData)
        FreqData(i).edge1 = FreqData(i).Value * 2^(0.5*log2((FreqData(i-1).Value)/FreqData(i).Value));
        FreqData(i).edge2 = FreqData(i).Value * 2^(-0.5*log2((FreqData(i-1).Value)/FreqData(i).Value)); %This assures the upper half of the first rectangle has the same size of the lower half (which has data to determine it)
    end
end

%% Main calculations (3) - Evaluating and summing the contribution of each frequency component to each band
%Preparing variables for GPU computing form
FreqDataEdge1 = zeros(1,length(FreqData));
FreqDataEdge2 = zeros(1,length(FreqData));
FreqDataValue = zeros(1,length(FreqData));

bandCenter = zeros(1,length(band));
bandValue = zeros(1,length(band));
bandLength = length(band);


for j = 1:length(FreqData)
    FreqDataEdge1(j) = FreqData(j).edge1;
    FreqDataEdge2(j) = FreqData(j).edge2;
    FreqDataValue(j) = FreqData(j).Value;
end
for i = 1:length(band)
    bandCenter(i) = band(i).center;
end
FreqDataLength = length(FreqData);


for i = 1:bandLength
    
    if i > 1 & i < bandLength 
        centers = [bandCenter(i-1) bandCenter(i) bandCenter(i+1)]; %Smaller band centers vector for enhanced computational speed
    end
    if i == 1
        centers = [0 bandCenter(i) bandCenter(i+1)]; %Smaller band centers vector for enhanced computational speed
    end
    if i == bandLength 
        centers = [bandCenter(i-1) bandCenter(i)]; %Smaller band centers vector for enhanced computational speed
    end
    
    for j = 1:FreqDataLength
        %Step1: Check which frequency inputs fall inside the band being evaluated
        %Step2: Correct boundaries f1 and f2 for calculation if necessary
        %Step3: Convert f1 and f2 to n values
        %Step4: Sum to the band the contribution of that frequency
        %       compenent, which will be Value*Coefficient - provided by 
        %       integrating the curve formula through the bandwidth 
        %       occupied by the frequency component in avg mode, and by the
        %       curve formula itself in sum mode.
        
        partialBandValue = zeros(1,1);
        
        if i > 1 & i < bandLength
            %S1
            if eq(or(ge(FreqDataEdge1(j),centers(3)),le(FreqDataEdge2(j),centers(1))),false)
                %S2
                f1 = FreqDataEdge1(j); %Initial values for f1 and f2
                f2 = FreqDataEdge2(j); %Initial values for f1 and f2
                if (f1 < centers(1))  %Correcting if necessary
                    f1 = centers(1);
                end
                if (f2 > centers(3)) %Correcting if necessary
                    f2 = centers(3);
                end
                
                %S3
                n_f1 = b*log2(f1/centers(2));
                n_f2 = b*log2(f2/centers(2));
                n_fc = b*log2(FreqDataValue(j)/centers(2));
                
                %S4
                    if Mode == 'avg'
                        coef = (0.5*n_f2 + sin(pi*n_f2)/(2*pi)) - (0.5*n_f1 + sin(pi*n_f1)/(2*pi)); % Determining fraction of contribution of current input frequency(j) to the band according to the model used
                    end
                    if Mode == 'sum'
                        coef = (cos(0.5*n_fc*pi))^2;
                    end
                    partialBandValue = ValueData(j)*coef; %Individual, partial, contribution to band value evaluated
            end
        end
       
        
        if i == 1
            %S1
            if eq(or(ge(FreqDataEdge1(j),centers(3)),le(FreqDataEdge2(j),evaluationedge1)),false) %Makes sure the components near the 1st band are inside the limits of the band
                %S2
                f1 = FreqDataEdge1(j); %Initial values for f1 and f2
                f2 = FreqDataEdge2(j); %Initial values for f1 and f2
                if (f1 < evaluationedge1)  %Correcting if necessary
                    str = strcat('Data was ignored at the very low frequency region! Cutoff frequency in use:',32,num2str(evaluationedge1),'Hz. For frequency',32,num2str(FreqDataValue(j)),'Hz, its window would start at',32, num2str(FreqDataEdge1(j)),'Hz.');
                    warning(str)
                    f1 = evaluationedge1;
                end
                if (f2 > centers(3)) %Correcting if necessary
                    f2 = centers(3);
                end
                
                %S3
                n_f1 = b*log2(f1/centers(2));
                n_f2 = b*log2(f2/centers(2));
                
                %S4
                coef = (0.5*n_f2 + sin(pi*n_f2)/(2*pi)) - (0.5*n_f1 + sin(pi*n_f1)/(2*pi)); % Determining fraction of contribution of current input frequency(j) to the band according to the model used
                partialBandValue = ValueData(j)*coef; %Individual, partial, contribution to band value evaluated
            end
        end
            
            
        if i == bandLength
            %S1
            if eq(or(ge(FreqDataEdge1(j),evaluationedge2),le(FreqDataEdge2(j),centers(1))),false) %Makes sure the components near the last band are inside the limits of the band
                %S2
                f1 = FreqDataEdge1(j); %Initial values for f1 and f2
                f2 = FreqDataEdge2(j); %Initial values for f1 and f2
                if (f1 < centers(1))  %Correcting if necessary
                    f1 = centers(1);
                end
                if (f2 > evaluationedge2) %Correcting if necessary
                    str = strcat('Data was ignored at the very high frequency region! Cutoff frequency in use:',32,num2str(evaluationedge2),'Hz. For frequency',32,num2str(FreqDataValue(j)),'Hz, its window would start at',32, num2str(FreqDataEdge2(j)),'Hz.');
                    warning(str)
                    f2 = evaluationedge2;
                end
                %S3
                n_f1 = b*log2(f1/centers(2));
                n_f2 = b*log2(f2/centers(2));
                
                %S4
                coef = (0.5*n_f2 + sin(pi*n_f2)/(2*pi)) - (0.5*n_f1 + sin(pi*n_f1)/(2*pi)); % Determining fraction of contribution of current input frequency(j) to the band according to the model used
                partialBandValue = ValueData(j)*coef; %Individual, partial, contribution to band value evaluated
            end
        end
        bandValue(i) = partialBandValue + bandValue(i); %Band value being updated
    end
end
% inputFreq
% ValueData
% bandCenter
% bandValue
% error('Debugging stop')

if plots == true
    figure(99)
    ax = gca;
    scatter(inputFreq,ValueData)
    hold on
    plot(bandCenter,bandValue,'-*')
    xlim([inputFreq(1) inputFreq(end)])
    xlim([20 20000])
    ylim([min(ValueData) max(ValueData)])
    ax.XScale = 'log';
    hold off
end

end


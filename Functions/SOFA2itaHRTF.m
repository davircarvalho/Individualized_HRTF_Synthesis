function result = SOFA2itaHRTF(handleSofa)

%% check HRTF
isHrtf = ismember(handleSofa.GLOBAL_SOFAConventions,{'SimpleFreeFieldHRIR','SimpleFreeFieldHRTF','SimpleFreeFieldTF'});
if ~isHrtf
    error("%s: Datatype not typically used for HRTFs",handleSofa.GLOBAL_SOFAConventions);
end

%% get positional data
% head position coordinates and orientation
recieverCoordinates = ita_sofa_getCoordinates(handleSofa,'channelCoordinateType','ReceiverPosition');
listenerView        = ita_sofa_getCoordinates(handleSofa,'channelCoordinateType','ListenerView');
listenerUp          = listenerView;
listenerUp.(listenerUp.coordSystem) = handleSofa.ListenerUp;

%source coordinates
sourceCoordinates = ita_sofa_getCoordinates(handleSofa,'channelCoordinateType','SourcePosition');

%% get datatype info
switch handleSofa.GLOBAL_DataType
    case 'TF'
        audioData = getSofaFreqData(handleSofa);
    case {'FIR','FIR-E','FIRE'}
        audioData = getSofaTimeData(handleSofa);
    case 'SOS'
        error('SOS type not yet implemented for ITA-Toolbox');
end

% set direction
audioData(1).channelCoordinates = sourceCoordinates;
audioData(2).channelCoordinates = sourceCoordinates;

%% final data as HRTF
if isa(audioData,'itaResult')
    result = audioData;
    return
else
    result = itaHRTF(audioData);
end

result.objectCoordinates = recieverCoordinates;
result.objectUpVector = listenerUp;
result.objectViewVector = listenerView;

%% add metaData as userData
result.userData = getSofaMetadata(handleSofa);
end

function audioData = getSofaTimeData(handleSofa)
data = handleSofa.Data.IR;
samplingrate = handleSofa.Data.SamplingRate;
if size(data,2) ~= 2
    error('unknown data structure')
end
leftEarData = itaAudio(squeeze(data(:,1,:))',samplingrate,'time');
rightEarData = itaAudio(squeeze(data(:,2,:))',samplingrate,'time');
audioData = [leftEarData, rightEarData];
end

function [audioData] = getSofaFreqData(handleSofa,sArgs)
% handle frequency domain transfer functions:
%   decide wether to return itaAudio or itaResult object
if ~ismember(handleSofa.N_LongName,{'f','frequency'})
    error('Unspecified content of N in sofa object');
end
if ~ismember(handleSofa.N_Units,{'hertz','hz','Hertz','Hz'})
    error('Uknown frequency unit in sofa object');
end

frequencies = handleSofa.N;
data = complex(handleSofa.Data.Real,handleSofa.Data.Imag);

if std(diff(frequencies)) > 0.01*mean(diff(frequencies))
    ita_verbose_info('Please check frequency Vector. It seems there are some values missing - returnin itaResult',0);
    dataL =  itaResult(squeeze(data(:,1,:)).',frequencies,'freq');
    dataR =  itaResult(squeeze(data(:,2,:)).',frequencies,'freq');
    audioData = [dataL, dataR];
    return
else
    ita_verbose_info('Automatically converting from equidistantly sampled frequency domain data, handle with care',1)
end


samplingrate = 2*max(frequencies);
% check number of ears
switch size(data,2)
    case 1
        ita_verbose_info('ita_read_sofa_hrtf: Only data of one ear present, using it for both ears',0)
        dataL =  itaAudio(squeeze(data(:,1,:)).',samplingrate,'freq');
        dataR =  dataL;
    case 2 
        dataL =  itaAudio(squeeze(data(:,1,:)).',samplingrate,'freq');
        dataR =  itaAudio(squeeze(data(:,2,:)).',samplingrate,'freq');
    otherwise
        error('Number of ears not matching')
end

audioData = [dataL, dataR];
end

function result = getSofaMetadata(handleSofa)
% get meta data to keep with the ITA files

%list of possible dataFields
metaDataFields = {'GLOBAL_Conventions','GLOBAL_Version','GLOBAL_SOFAConventions','GLOBAL_SOFAConventionsVersion' ...
    ,'GLOBAL_APIName','GLOBAL_APIVersion','GLOBAL_ApplicationName','GLOBAL_ApplicationVersion','GLOBAL_AuthorContact' ...
    ,'GLOBAL_Comment','GLOBAL_DataType','GLOBAL_History','GLOBAL_License','GLOBAL_Organization','GLOBAL_References' ...
    ,'GLOBAL_RoomType','GLOBAL_Origin','GLOBAL_DateCreated','GLOBAL_DateModified','GLOBAL_Title','GLOBAL_DatabaseName' ...
    ,'GLOBAL_RoomDescription','GLOBAL_ListenerShortName','API','ListenerPosition','ListenerPosition_Type','ListenerPosition_Units'...
    ,'EmitterPosition','EmitterPosition_Type','EmitterPosition_Units','RoomCornerA','RoomCornerA_Type','RoomCornerA_Units' ...
    ,'RoomCornerB','RoomCornerB_Type','RoomCornerB_Units','','','','','','',''};

for index = 1:length(metaDataFields)
    if isfield(handleSofa,metaDataFields{index})
        result.(metaDataFields{index}) =  handleSofa.(metaDataFields{index});
    end
end

end %getSofaMetadata()
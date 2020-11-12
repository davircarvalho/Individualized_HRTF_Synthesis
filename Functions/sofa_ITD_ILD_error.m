function [ITD_error, ILD_error] = sofa_ITD_ILD_error(Obj_med, Obj_ref, varargin)
% optional parameters:
%         'time':       ITD error output will have units in seconds
%         'samples'     ITD error output will have units in samples

%% Parse inputs
% Opção de saída em samples ou tempo
defaultMode = 'time';
validOutputs = {'samples','time'};
checkOutMode = @(x) any(validatestring(x, validOutputs));

p = inputParser;
addRequired(p,'Obj_med',@isstruct);
addRequired(p,'Obj_ref',@isstruct);
addOptional(p,'outputMode', defaultMode,checkOutMode)
parse(p, Obj_med, Obj_ref, varargin{:})


%%
    IR_ref = shiftdim(Obj_ref.Data.IR, 2);
    IR_med = shiftdim(Obj_med.Data.IR, 2);
    N = size(IR_ref, 1);
    
%% ILD --------------------------------------------------------------------
    % diferenca de fase interaural
    ITF_ref  = fft(IR_ref(:,:,1), N)./fft(IR_ref(:,:,2), N);
    ITF_msrd = fft(IR_med(:,:,1), N)./fft(IR_med(:,:,2), N);
    % ild
    ILD_ref  = 20*log10(abs(ITF_ref(2:N/2,:)));
    ILD_msrd = 20*log10(abs(ITF_msrd(2:N/2,:)));
    % error
    ILD_error = mean(abs(ILD_ref - ILD_msrd));
    
%% ITD --------------------------------------------------------------------
    ITD_ref = sofaGetITD(Obj_ref, p.Results.outputMode);
    ITD_med = sofaGetITD(Obj_med, p.Results.outputMode);    
    % error
    ITD_error = abs(ITD_ref - ITD_med);

%     figure()
%     plot(ITD_ref); hold on
%     plot(ITD_med); 
%     hold off
%     legend('ref', 'sim')
end 
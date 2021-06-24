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
    
%% ILD --------------------------------------------------------------------
    IR_ref = shiftdim(Obj_ref.Data.IR, 2);
    IR_med = shiftdim(Obj_med.Data.IR, 2);
    
    ILD_ref = getILD(IR_ref(:,:,1), IR_ref(:,:,2));
    ILD_med = getILD(IR_med(:,:,1), IR_med(:,:,2));

    % error
    ILD_error = abs(ILD_ref - ILD_med);
    
%% ITD --------------------------------------------------------------------
    ITD_ref = SOFAgetITD(Obj_ref, p.Results.outputMode, 'thr', 20);
    ITD_med = SOFAgetITD(Obj_med, p.Results.outputMode, 'thr', 20); 
    
    % error
    ITD_error = abs(ITD_ref - ITD_med);

%     figure()
%     plot(ITD_ref); hold on
%     plot(ITD_med); 
%     hold off
%     legend('ref', 'sim')
end 


function ILD = getILD(L, R)
pLeft  = fft(L);
pRight = fft(R);

pLeft(1,:) = [];
pRight(1,:) =[];

pLeftIntegral  = sum(abs(pLeft).^2, 1 );
pRightIntegral = sum(abs(pRight).^2, 1 );

ILD = 10*log10( pLeftIntegral ./ pRightIntegral );
end






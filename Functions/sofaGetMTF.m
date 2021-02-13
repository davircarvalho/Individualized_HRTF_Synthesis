function mtf = sofaGetMTF(Obj)
% Determinação do ITD para HRTFs SOFA
% Davi R. Carvalho @UFSM - Engenharia Acustica - Setembro/2020
    %% Comum 
    mtf = zeros(length(Obj.SourcePosition), 1);
    IR = shiftdim(Obj.Data.IR, 2);
    fs = Obj.Data.SamplingRate;
    azim = 0; elev = 0;
    pos = Obj.SourcePosition;
    [~,idx_pos] = min(sqrt((pos(:,1)-azim).^2 + (pos(:,2)-elev).^2));
    
    %% Metodo 1 : Threshold (usando ita)
    % a partir da diferença entre picos de RI (ISO 3382 A.3.4.).
    A = itaAudio; A.samplingRate = fs;
    B=A;
    A.time = IR(:,idx_pos,1); % L
    for k = 1:size(IR, 2)      
        B.time = IR(:,k,1); % L    
         % ponto de chegada 
        OnSetL = ita_start_IR(A, 'correlation', false, 'threshold', 10);
        OnSetR = ita_start_IR(B, 'correlation', false, 'threshold', 10);
        
         % diferenca interaural
        mtf(k) = abs(OnSetL - OnSetR);
    end 
end
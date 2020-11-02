function Obj = sofaNormalize(Obj)
% Normalizes channel balance acording to mismatch between level for the 
% source right in front of the subject and checks if the left channel is 1 
% and right channel is 2, if misplaced it'll be corrected
elev = 0; azim = 0;
posi = round(Obj.SourcePosition);
idx_pos = dsearchn(posi(:,[1,2]),  [azim, elev]);

% normalize balance
norm =  max(abs(Obj.Data.IR(idx_pos, :,:)),[],3);
IR(:,1,:) = Obj.Data.IR(:,1,:)./norm(1);
IR(:,2,:) = Obj.Data.IR(:,2,:)./norm(2);

% normalize peaks
IR = IR./max(abs(IR(:)));

%% Check IR channels
azim = 90;
idx_pos = dsearchn(posi(:,[1,2]),  [azim, elev]);
chkIR =  max(abs(Obj.Data.IR(idx_pos, :,:)),[],3);
if chkIR(1)  > chkIR(2) % espera-se que ipsislateral sempre tenha maior amplitude
    Obj.Data.IR = IR;
else % exchange posi
    % double check
    idx_pos = dsearchn(posi(:,[1,2]),  [270, elev]);
    chkIR =  max(abs(Obj.Data.IR(idx_pos, :,:)),[],3);
    if chkIR(1)  > chkIR(2) 
        warning('on')
        warning('left and right IRs are now channel 1 and 2, as SOFA convention 1.0 expects')
        Obj.Data.IR(:,1,:) = IR(:,2,:);
        Obj.Data.IR(:,2,:) = IR(:,1,:);
    end
end


end
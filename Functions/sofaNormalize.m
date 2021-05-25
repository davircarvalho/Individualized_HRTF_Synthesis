function Obj = sofaNormalize(Obj)
% Normalizes channel balance acording to mismatch between level for the 
% source right in front of the subject 
elev = 0; azim = 0;
posi = round(Obj.SourcePosition);
idx_pos = dsearchn(posi(:,[1,2]),  [azim, elev]);

% normalize L/R balance
norm =  max(abs(Obj.Data.IR(idx_pos, :,:)),[],3);
IR(:,1,:) = Obj.Data.IR(:,1,:)./norm(1);
IR(:,2,:) = Obj.Data.IR(:,2,:)./norm(2);

% normalize peaks
IR = IR./max(abs(IR(:)));


%% Check IR channels (deprecated)
% This was build before the correction in the CIPIC to SOFA fix
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
        warning(['Left and right channels were inverted in order to '...
                 'respect ipsislateral/contralateral expected behaviour. ' ...
                 'Be aware that the coordinates could be the real source of '...
                 'the mismatch, but that cannot be identified by the sole evaluation, '...
                 'of the SOFA object.'])
        Obj.Data.IR(:,1,:) = IR(:,2,:);
        Obj.Data.IR(:,2,:) = IR(:,1,:);
    end
end


end
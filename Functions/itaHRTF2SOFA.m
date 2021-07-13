function Obj = itaHRTF2SOFA(this)
    Obj = SOFAgetConventions('SimpleFreeFieldHRIR');
    
    % IR
    Obj.Data.SamplingRate = this.samplingRate;
    irL = this.getEar('L');
    irR = this.getEar('R');
    Obj.Data.IR = nan(this.dirCoord.nPoints, 2, size(irL.time, 1));
    Obj.Data.IR(:,1,:) = (irL.time).';
    Obj.Data.IR(:,2,:) = (irR.time).';
    
    % Coordinates
    pos(:,1) = this.dirCoord.phi_deg;
    pos(:,2) = this.dirCoord.theta_deg -90;
    pos(:,3) = this.dirCoord.r;
    Obj.SourcePosition = pos;
    Obj = SOFAupdateDimensions(Obj);
end
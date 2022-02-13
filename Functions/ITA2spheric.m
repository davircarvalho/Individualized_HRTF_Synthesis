function Obj = ITA2spheric(Obj)
    for l = 1:length(Obj.SourcePosition)
        x = Obj.SourcePosition(l, 1);  
        y = Obj.SourcePosition(l, 2); 
        z = Obj.SourcePosition(l, 3);
        % new coordinates
        [az,elev,r] = cart2sph(x,y,z);
        azi=rad2deg(az); elev=rad2deg(elev);
        [azi,ele]   = nav2sph(azi,elev);
        azi(azi == 360) = 0;
        % update coordinates
        Obj.SourcePosition(l, 1) = azi;
        Obj.SourcePosition(l, 2) = ele; 
        Obj.SourcePosition(l, 3) = round(r);
        % more metadata
        Obj.SourcePosition_Type = 'spherical';
        Obj.SourcePosition_Units = 'degree, degree, meter';              
    end  
end
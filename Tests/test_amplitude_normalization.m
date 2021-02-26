% 
clear all; clc
% find distance between th source position and each ear 
Obj = SOFAload('individuo_141.sofa');

%
L_ear = Obj.ReceiverPosition(1,:);
R_ear = Obj.ReceiverPosition(2,:);

% 
pos = Obj.SourcePosition(100,:);
[tx, ty, tz] = sph2cart(pos(1), pos(2), pos(3));

%
Ldist = sqrt((L_ear(1) - tx)^2 + (L_ear(2) - ty)^2 + (L_ear(3) - tz)^2);
Rdist = sqrt((R_ear(1) - tx)^2 + (R_ear(2) - ty)^2 + (R_ear(3) - tz)^2);

% Calculate distance normalization
dist_norm = ((0.5/Ldist) + (0.5/Rdist)).^-1;
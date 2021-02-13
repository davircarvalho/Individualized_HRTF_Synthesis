%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function [samplingGridParL, samplingGridParR] = supdeq_parallax(samplingGrid, sourceDistance, radius, earPosition) 
%
% This function applies parallax shifts to a (far-field) sampling grid according 
% to the source distance and the ear position, and returns the new sampling 
% grid including parallax shifts for the left and right ear.
%
% Output:
% samplingGridParL/R    - Parallax-shifted sampling grid for left and right ear
%
% Input:
% samplingGrid          - Spatial sampling grid (Q x 2 matrix), where the first 
%                       column holds the azimuth and the second the elevation (both in degree).
%                       Azimuth in degree (0=front, 90=left, 180=back, 270=right)
%                       (0 points to positive x-axis, 90 to positive y-axis)
%                       Elevations in degree (0=North Pole, 90=front, 180=South Pole)
%                       (0 points to positive z-axis, 180 to negative z-axis)
% sourceDistance        - Distance (in meter) to sound source
% radius                - Head/Sphere radius in m
% earPosition           - 4 x 1 row vector describing the position of the ears in
%                       spherical coordinates in degree [azL, elL, azR, elR]
%                       Default: [90, 90, 270, 90] (left-right symmetrical)
%
% Dependencies: -
%
% (C) 2020 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function [samplingGridParL, samplingGridParR] = supdeq_parallax(samplingGrid, sourceDistance, radius, earPosition)

if nargin < 4 || isempty(earPosition)
    earPosition = [90, 90, 270, 90];
end

if sourceDistance < radius
    error('Source distance smaller than head radius!');
end

%% Transform sampling grid and ear positions to radiant in matlab coordinates system

samplingGridRad(:,1:2) = samplingGrid(:,1:2)*pi/180;
samplingGridRad(:,2) = pi/2 - samplingGridRad(:,2);

earPosition_L = earPosition(1:2)*pi/180;
earPosition_L(2) = pi/2 - earPosition_L(2);

earPosition_R = earPosition(3:4)*pi/180;
earPosition_R(2) = pi/2 - earPosition_R(2);

%% Transform sampling grid and ear positions to cartesian coordinates and apply parallax shift

[x,y,z] = sph2cart(samplingGridRad(:,1),samplingGridRad(:,2),sourceDistance);

[x_ear_L,y_ear_L,z_ear_L] = sph2cart(earPosition_L(:,1),earPosition_L(:,2),radius); 
[x_ear_R,y_ear_R,z_ear_R] = sph2cart(earPosition_R(:,1),earPosition_R(:,2),radius);  

xL = x-x_ear_L;
xR = x-x_ear_R;
yL = y-y_ear_L;
yR = y-y_ear_R;
zL = z-z_ear_L;
zR = z-z_ear_R;

[AzParL,ElParL] = cart2sph(xL,yL,zL);
[AzParR,ElParR] = cart2sph(xR,yR,zR);

%% Transform back to degree in SH coordinates system

AzParL = mod(AzParL,2*pi);
AzParR = mod(AzParR,2*pi);
AzParL = AzParL*180/pi;
AzParR = AzParR*180/pi;

ElParL=pi/2-ElParL;
ElParR=pi/2-ElParR;
ElParL = ElParL*180/pi;
ElParR = ElParR*180/pi;

%% Store in new samplingGrids for left and right channel

samplingGridParL = [AzParL,ElParL];
samplingGridParR = [AzParR,ElParR];

end


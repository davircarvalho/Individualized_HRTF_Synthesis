function [y,z1] = firfilter(b,x)
% Inputs:
%   b - 1xNTaps row vector of coefficients
%   x - A frame of  noisy input 

% States:
%   z, z1 - NTapsx1 column vector of states

% Output:
%   y - A frame of filtered output
 
persistent z;

if (isempty(z))
    z = zeros(length(b),1);
end
Lx = size(x,1);
y = zeros(size(x),'like',x);

z1 = z;
for m = 1:Lx
    % Load next input sample
    z1(1,:) = x(m,:);
    
    % Compute output
    y(m,:) = b*z1;
    
    % Update states
    z1(2:end,:) = z1(1:end-1,:);
    z = z1;
end
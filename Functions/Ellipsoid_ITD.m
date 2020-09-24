function itd = Ellipsoid_ITD(wid, dep, hei, pos)
% inputs idw = width 
%        dep = depth
%        hei = height
%        pos = positions
%% 

for n = 1:length(pos)  
    r_ellipseLL(n) = wid./sqrt(1 - (((d_ellipse^2 - wid^2)/d_ellipse^2)*cos(pos(n,1))).^2);
end
        
 
freq = 20:20000; % frequencias
c0 = 343; % velocidade do som 
k = 2*pi*freq/c0;
a = mean(r_ellipse);
N_order = 10; % ordem do polinomio de legendre

for n = 1:length(pos)
    for q = 1:N_order
      temp(q) = 1j^(q+1)*legendreP(q,cos(n))*(2*q+1)./(besselj(q, k*a) - 1j* bessely(k*a));    
    end
    temp2(n) = sum(temp);
end
TFspere_l = 1/(k.*a).^2 * temp(n);



end
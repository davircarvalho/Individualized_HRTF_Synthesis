function [IR1, IR2] = hrtf_simetry_split(Obj)
% EXPLORAR REDUNDANCIA DE DAS HRTFs
% Considerando que para angulos simetricos em relação a 0°
% (ex.: -30° e 30°)
% a resposta observada para cada orelha será virtualmente a mesma.
%       (ex.: a resposta vista pela orelha esquerda a um angulo de -30°
%       terá as mesmas propriedades da vista pela orelha direita a um angulo de 30°).
%


%%% INPUT %%%
% SOFA object  (PS.: make sure the coordinate system is spherical
%                                                azi: 0° -> 360°
%                                                ele: -90° -> 90°). 
%% 
posi1 = (Obj.SourcePosition);
IR1   = Obj.Data.IR(:,1,:);

%% Find corresponding positions
idx = zeros(length(posi1), 1);
for k = 1:length(posi1)
    if posi1(k,1) ~= 0
        idx(k) = find(posi1(:,1) == (360-posi1(k,1)) & posi1(:,2) == posi1(k,2));
    else % azimute = 0°
        idx(k) = k;
    end 
end    

% Address new positions
posi2 = Obj.SourcePosition(idx, :);
[posi2, idx_sort] = sortrows( posi2, [-2,1]);

IR2 = Obj.Data.IR(idx_sort, 2, :);

%% PREPARE OUTPUT and metadata
clc
Obj1 = SOFAgetConventions('SimpleFreeFieldHRIR');
Obj1.Data.IR = IR1; 
Obj1.SourcePosition = posi1;
Obj1.Data.SamplingRate = Obj.Data.SamplingRate;


Obj2 = SOFAgetConventions('SimpleFreeFieldHRIR');
Obj2.Data.IR = IR2; 
Obj2.SourcePosition = posi2;
Obj2.Data.SamplingRate = Obj.Data.SamplingRate;

Obj1.API.R = 1;
Obj2.API.R = 1;

%% Plot symmetri between channels

% type = 'MagHorizontal';
% subplot(211)
% SOFAplotHRTF(Obj1, type);
% subplot(212)
% SOFAplotHRTF(Obj2, type);

%% PLOT - HRIR [Simulada vs Medida]
% figure()
% 
% % DEFINA A DIREÇÃO DA RI
% azim = 20; 
% elev = 10;
% 
% % Get index of measurements with the same directions
% pos1=find(round(Obj1.SourcePosition(:,1))==azim & round(Obj1.SourcePosition(:,2))==elev);
% pos2=find(round(Obj2.SourcePosition(:,1))==azim & round(Obj2.SourcePosition(:,2))==elev);
% 
% %%% ESQUERDA --------------------------------------------------------------
% % Medida
% plot(real(squeeze(Obj1.Data.IR(pos1, 1, :))), 'linewidth', 1.5);hold on
% % Proposta
% plot(real(squeeze(Obj2.Data.IR(pos2, 1, :))), 'r', 'linewidth', 1.5); hold off
% 
% out_pos = Obj.SourcePosition;
% subject = 1;
% xlabel('Amostras')
% ylabel('Esquerda')
% azim = num2str(out_pos(pos1,1));
% elee = num2str(out_pos(pos1,2));
% title(['HRIR Invididuo ', num2str(subject), ', azimute ', azim, ...
%                                    '°, elevação ', elee, '°']);
% legend('L', 'R')
% set(gca,'FontSize',12)
% axis tight

end
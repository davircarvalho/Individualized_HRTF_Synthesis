function [hrir_l,hrir_r,name]=AA_ARI2CIPIC(hM,meta,stimPar,showfigures)
% AA_ARI2CIPIC - Convert ARI HRTF data to CIPIC HRTF format
% 
% [hrir_l,hrir_r,name]=AA_ARI2CIPIC(hM,meta,stimPar,showfigures)
%
% Input:
%  hM: ARI HRTF data
%  meta: hM structured meta data
%  stimPar: stimulation parameters
%  showfigures: 0: no figures, ~=0 show figures
% Output:
%  hrir_l: CIPIC HRTF data for left ear
%  hrir_r: CIPIC HRTF data for right ear
%  name: Subject ID
% 
% (see 'ARI HRTF format doc' and 'Documentation for the UCD HRIR Files' for further details)
% 
% by Michael Mihocic 12.05.2011
%  Austrian Academy of Sciences, Acoustics Research Institute
% Last change: 05.10.2011 by Michael Mihocic, upgrade to structured meta
%  data hM version 2.0.0

if isempty(showfigures)
    showfigures=1;
end

% CIPIC original positions (optimal values)
azi=[-80 -65 -55 -45:5:45 55 65 80]; % CIPIC Original Azi
ele=-45+5.625*(0:49); % CIPIC Ele

pos=nan(length(azi)*length(ele),12); 
% pos columns
%  Col 1: CIPIC Lat [-90 to 90], horizontal-polar
%  Col 2: CIPIC Pol [-90 to 180], horizontal-polar
%  Col 3: CIPIC Azi [0 to 359], geodetic
%  Col 4: CIPIC Ele [-90 to 90], geodetic
%  Col 5: CIPIC Lat rounded to ARI Grid, horizontal-polar
%  Col 6: CIPIC Pol rounded to ARI Grid, horizontal-polar
%  Col 7: CIPIC Azi rounded to ARI Grid, geodetic
%  Col 8: CIPIC Ele rounded to ARI Grid, geodetic
%  Col 9: Lat Error (CIPIC Original / rounded in ARI Grid), horizontal-polar
%  Col 10: Pol Error (CIPIC Original / rounded in ARI Grid), horizontal-polar
%  Col 11: Azi Error (CIPIC Original / rounded in ARI Grid), geodetic
%  Col 12: Ele Error (CIPIC Original / rounded in ARI Grid), geodetic

idx=zeros(length(azi)*length(ele),1);
zz=1;

% Get list of all CIPIC positions
for ii=1:length(azi)    
    for jj=1:length(ele)
        pos(zz,1)=azi(ii);
        pos(zz,2)=ele(jj);    
        
%         Find nearest position of CIPIC data in ARI grid (meta.pos from input)
        tsqr=sqrt((meta.pos(:,6)-pos(zz,1)).^2+(meta.pos(:,7)-pos(zz,2)).^2);
        idx(zz)=find(tsqr==min(tsqr),1);
        [pos(zz,3),pos(zz,4)]=hor2geo(pos(zz,1),pos(zz,2));
        pos(zz,5)=meta.pos(idx(zz),6);
        pos(zz,6)=meta.pos(idx(zz),7);
        [pos(zz,7),pos(zz,8)]=hor2geo(pos(zz,5),pos(zz,6));
        zz=zz+1;
    end
end

% if showfigures ~= 0 
%     % Plot Comparison CIPIC Original and rounded to ARI grid
%     
%     % Error values only needed of plotted
%     pos(:,9)=abs(abs(pos(:,1))-abs(pos(:,5)));
%     pos(:,10)=abs(abs(pos(:,2))-abs(pos(:,6)));
%     pos(:,11)=abs(abs(pos(:,3))-abs(pos(:,7)));
%     pos(:,12)=abs(abs(pos(:,4))-abs(pos(:,8)));
%     
%     % sort...
%     uni_lat=unique(pos(:,1));
%     uni_pol=unique(pos(:,2));
% 
%     % define label values
%     azilabel=[-80 -40 -20 0 20 40 80];
%     elelabel=[-45 0 45 90 135 180 225];
% 
%     kk=1; ll=1;
%     azitick=zeros(size(azilabel));
%     eletick=zeros(size(elelabel));
% 
%     err_lat=nan(size(uni_pol,1),size(uni_lat,1));
%     err_pol=nan(size(uni_pol,1),size(uni_lat,1));
%     err_azi=nan(size(uni_pol,1),size(uni_lat,1));
%     err_ele=nan(size(uni_pol,1),size(uni_lat,1));   
%     
% %     Calculate errors
%     for ii=1:size(uni_lat,1)
%         for jj=1:size(uni_pol,1)
%             if ~isempty(find(pos(:,1)==uni_lat(ii) & pos(:,2)==uni_pol(jj), 1));
%                 err_lat(jj,ii)=pos(pos(:,1)==uni_lat(ii) & pos(:,2)==uni_pol(jj) ,9);
%                 err_pol(jj,ii)=pos(pos(:,1)==uni_lat(ii) & pos(:,2)==uni_pol(jj) ,10);
%                 err_azi(jj,ii)=pos(pos(:,1)==uni_lat(ii) & pos(:,2)==uni_pol(jj) ,11);
%                 err_ele(jj,ii)=pos(pos(:,1)==uni_lat(ii) & pos(:,2)==uni_pol(jj) ,12);                
%             end
%             if ii==1 % do only once
%                 if ~isempty(find(elelabel==uni_pol(jj,1), 1))
%                     eletick(ll)=jj;
%                     ll=ll+1;
%                 end
%             end
%         end
% 
%         if ~isempty(find(azilabel==uni_lat(ii,1), 1))
%             azitick(kk)=ii;
%             kk=kk+1;
%         end
%     end
% 
%     figure('Position',[1 30 600 500]);
%     
%     visdata([pos(:,7) pos(:,8)]);      % CIPIC in ARI Grid
%     visdata([pos(:,3) pos(:,4)],'rx'); % CIPIC Original
%     title 'CIPIC original (red) and ARI2CIPIC (blue) positions';
%     % % Save figures
%     % saveas(gcf, 'CIPIC ARI2CIPIC 0,90', 'fig');
%     % saveas(gcf, 'CIPIC ARI2CIPIC 0,90', 'bmp');
%     % view(0,0);
%     % saveas(gcf, 'CIPIC ARI2CIPIC 0,0', 'fig');
%     % saveas(gcf, 'CIPIC ARI2CIPIC 0,0', 'bmp');
%     view(45,45);
%     % saveas(gcf, 'CIPIC ARI2CIPIC 45,45', 'fig');
%     % saveas(gcf, 'CIPIC ARI2CIPIC 45,45', 'bmp');
%     % view(90,0);
%     % saveas(gcf, 'CIPIC ARI2CIPIC 90,0', 'fig');
%     % saveas(gcf, 'CIPIC ARI2CIPIC 90,0', 'bmp');
% 
%     % Plot Lateral Error
%     figure('Position',[8 531 560 420]);
%     pcolor(err_lat);
%     colorbar;
%     % shading flat;
%     ylabel('Pol CIPIC');
%     xlabel('Lat CIPIC');
%     set(gca,'XTick',azitick);
%     set(gca,'YTick',eletick);
%     set(gca,'YTickLabel',elelabel)
%     set(gca,'XTickLabel',azilabel)
%     title('Lateral Error [deg]');
% 
%     % Plot Polar Error
%     figure('Position',[588 531 560 420]);
%     pcolor(err_pol);
%     colorbar;
%     % shading flat;
%     ylabel('Pol CIPIC');
%     xlabel('Lat CIPIC');
%     set(gca,'XTick',azitick);
%     set(gca,'YTick',eletick);
%     set(gca,'YTickLabel',elelabel)
%     set(gca,'XTickLabel',azilabel)
%     title('Polar Error [deg]');
%     
%         % Plot Azimuth Error
%     figure('Position',[8 531 560 420]);
%     pcolor(err_azi);
%     colorbar;
%     % shading flat;
%     ylabel('Ele CIPIC');
%     xlabel('Azi CIPIC');
%     set(gca,'XTick',azitick);
%     set(gca,'YTick',eletick);
%     set(gca,'YTickLabel',elelabel)
%     set(gca,'XTickLabel',azilabel)
%     title('Azimuth Error [deg]');
% 
%     % Plot Elevation Error
%     figure('Position',[588 531 560 420]);
%     pcolor(err_ele);
%     colorbar;
%     % shading flat;
%     ylabel('Ele CIPIC');
%     xlabel('Azi CIPIC');
%     set(gca,'XTick',azitick);
%     set(gca,'YTick',eletick);
%     set(gca,'YTickLabel',elelabel)
%     set(gca,'XTickLabel',azilabel)
%     title('Elevation Error [deg]');
% end

% Building hrir_l and hrir_r
hrir_l=zeros(length(azi),length(ele),length(hM(:,1,1)));
hrir_r=hrir_l;
for ii=1:length(pos)
    idx=find(meta.pos(:,6)==pos(ii,5) & meta.pos(:,7)==pos(ii,6));
    hrir_l(azi==pos(ii,1),ele==pos(ii,2),:)=hM(:,idx(1),1);
    hrir_r(azi==pos(ii,1),ele==pos(ii,2),:)=hM(:,idx(1),2);
end

% name
if isfield(stimPar,'SubjectID')
    name=stimPar.SubjectID;
else
    name='';
end
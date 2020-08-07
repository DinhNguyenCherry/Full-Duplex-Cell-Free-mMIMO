function [ done ] = Plot_Layout( RadiusOfRegion, Objects, Styles, Legends, WithColorBar )
%PLOT_LAYOUT Summary of this function goes here
%   Detailed explanation goes here

if (nargin<5)
    NoObj = length(Objects);
    WithColorBar = cell(1,NoObj);
    for iObj = 1:1:NoObj
        WithColorBar{1,iObj} = cell(1,1);
        WithColorBar{1,iObj}{1,1} = 'none';
    end
end



Extra = 0.05;

hold on

ang=0:0.01:2*pi; 

xp=RadiusOfRegion*cos(ang);
yp=RadiusOfRegion*sin(ang);
ax = plot(0+xp,0+yp, 'HandleVisibility', 'off');

% xp1=RadiusOfInnerzone*cos(ang);
% yp1=RadiusOfInnerzone*sin(ang);
% ax1 = plot(0+xp1,0+yp1);

NoObj = length(Objects);

for iObj = 1:1:NoObj
    
    if (strfind(WithColorBar{1,iObj}{1,1},'none'))
        scatter(Objects{1,iObj}(:,1), Objects{1,iObj}(:,2), 100, Styles{iObj})
    else
        scatter(Objects{1,iObj}(:,1), Objects{1,iObj}(:,2), 200, Objects{1,iObj}(:,3), 'filled')
        colorbar('Ticks',WithColorBar{1,iObj}{1,2},... 
                 'TickLabels',WithColorBar{1,iObj}{1,3})
             colormap(jet)
%         colorbar('Ticks',[min(Objects{1,iObj}(:,3)) 1 max(Objects{1,iObj}(:,3))],...
%                  'TickLabels',{'Sleep','Weak Active', 'Strong Active'})
        
%         AxesH = axes('CLim', [-12, 12]);
%         colorbar('YTick',0:1:1)
%         colorbar('YTickLabel', {'Sleep','Active'})
    end
%     scatter(positionUplinkUsers(:,1), positionUplinkUsers(:,2), 'b+ ')
%     scatter(positionDownlinkUsers(:,1), positionDownlinkUsers(:,2), 'r* ')
    
end


set(gca,'XTick',[-RadiusOfRegion :0.2*RadiusOfRegion: RadiusOfRegion])
set(gca,'YTick',[-RadiusOfRegion :0.2*RadiusOfRegion: RadiusOfRegion])

xlim([-(RadiusOfRegion+Extra*RadiusOfRegion) (RadiusOfRegion+Extra*RadiusOfRegion)])
ylim([-(RadiusOfRegion+Extra*RadiusOfRegion) (RadiusOfRegion+Extra*RadiusOfRegion)])
axis('square')

legend(Legends);

done = 1;


end


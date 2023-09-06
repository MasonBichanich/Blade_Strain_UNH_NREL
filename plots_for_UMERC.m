figure;
for i = 1:8
    h(i) = subplot_tight(8,1,i,[.05]);
    if i<5
        plot(vlink_8145{:,1},vlink_28175{:,i+1})
        ylim([-70 50])
        xlim([datetime('2022-11-22 19:45:18.156') datetime('2022-11-22 19:48:31.656')])
        %legend("SG - " + num2str(i))
        grid on
    else
        plot(vlink_28175{:,1},vlink_8145{:,i-3})
        ylim([-70 50])
        xlim([datetime('2022-11-22 19:45:18.156') datetime('2022-11-22 19:48:31.656')])
        %legend("SG - " + num2str(i))
        grid on
    end
end
hold on
set(gcf,'Color','w');
han=axes('position',[0.07 .275 .425 .525],'visible','off'); 
% han.Title.Visible='on';
% han.XLabel.Visible='on';
% han.YLabel.Visible='on';
ylabel('Microstrain','visible','on');
% xlabel(han,'Time [UTC]');
%% 
figure
polarplot(pi/180*phase_averaged{:,3}(:,1),(phase_averaged{:,3}(:,2)))
hold on
polarplot(0:.01:2*pi,zeros(length(0:.01:2*pi)),'k')
rlim([-50 50])
%%
figure
edges = pi/180*phase_averaged{:,3}(:,1);
counts = abs(phase_averaged{:,3}(:,2));
polarhistogram(BinEdges=[edges', 2*pi],BinCounts=counts')
rlim([-50 50])
%polarfill(gca,linspace(0,2*pi,length(phase_averaged{:,3}(:,2))),zeros(length(phase_averaged{:,3}(:,2))),phase_averaged{:,3}(:,2)-zeros(length(0:.01:2*pi)),phase_averaged{:,3}(:,2)+zeros(length(0:.01:2*pi)),'k',0.6)
%% ADV into different ensemble averages. This is not a real ensemble. I am cheating

vel = movmean(ADV2_xvel,64,'omitnan');
figure
good = ADV2_xvel < -1 & ADV2_xvel > -3;
plot(ADV2_datetime_UTC,vel)
xlim([datetime('2022-11-22 19:47:18.156') datetime('2022-11-22 19:48:15.656')])
yyaxis right
ax = gca;
ax.YColor = 'k';
plot(vlink_8145{:,1},vlink_8145{:,4},'k')
grid on
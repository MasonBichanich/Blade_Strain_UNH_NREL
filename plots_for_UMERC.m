
ax = mySubplot(8,1);
ylim([-70 50])
for i = 1:8
    if i<5
        plot(ax(i),vlink_8145{:,1},vlink_28175{:,i+1})
        set(ax(i),'ylim',[-70 50], ...
            'XTickLabel', '', ...
            'xlim', [datetime('2022-11-22 19:45:18.156') datetime('2022-11-22 19:48:31.656')]);
        grid(ax(i),'on')

    else
        plot(ax(i),vlink_28175{:,1},vlink_8145{:,i-3})
         set(ax(i),'ylim',[-70 50], ...
             'xlim', [datetime('2022-11-22 19:45:18.156') datetime('2022-11-22 19:48:31.656')]);
        grid(ax(i),'on')
         if i ~= 8
         set(ax(i), 'XTickLabel', '')
         end
    end 
end
hold on
set(gcf,'Color','w');
han=axes('position',[0.07 .275 .425 .525],'visible','off'); 
% han.Title.Visible='on';
% han.XLabel.Visible='on';
% han.YLabel.Visible='on';
ylabel('Microstrain','visible','on');
fontsize(16,'points')
linkaxes([ax(1) ax(2) ax(3) ax(4) ax(5) ax(6) ax(7) ax(8)],'x')
% xlabel(han,'Time [UTC]');
%%
ax = mySubplot(8,1);
plot(ax(1),1,1,'.')
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
ADV2_table = load('../../../Data_Files/MODAQ/Vector_Velocity_Data_2/Processed_Data/2022_11_22_to_11_22_Vector_2_UTC.mat');
ADV2_datenum_UTC  = table2array(ADV2_table.Vector_new(:,1));

ADV2_datetime_UTC = datetime(ADV2_datenum_UTC,'ConvertFrom','datenum','Format', 'MM/dd/yyyy HH:mm:ss.SSSSSSSSS');
ADV2_xvel = table2array(ADV2_table.Vector_new(:,2));
    vlink_8145 = readtable('../../../Data_Files/MODAQ/Vlink/Processed_Data/8145_vLink_strain-2022-11-22.csv');
    vlink_28175 = readtable('../../../Data_Files/MODAQ/Vlink/Processed_Data/28175_vLink_strain-2022-11-22.csv');

    upper_bound = 335729;  % this is where the turbine stops rotating
    % Extract vlink 8145 data
    vlink_8145 = vlink_8145(2:upper_bound,:);
    vlink_s1_s4_time = table2array(vlink_8145(:,1));
    vlink_s1_datenum = datenum(vlink_s1_s4_time);
    vlink_s1 = table2array(vlink_8145(:,2));
    vlink_s2 = table2array(vlink_8145(:,3));
    vlink_s3 = table2array(vlink_8145(:,4));
    vlink_s4 = table2array(vlink_8145(:,5));

    % Extract vlink 28175 data
    vlink_28175(1:2,:) = [];
    vlink_28175 = vlink_28175(1:upper_bound-1,:);
    vlink_s5_s8_time = table2array(vlink_28175(:,1));
    vlink_s5_datenum = datenum(vlink_s5_s8_time);
    vlink_s5 = table2array(vlink_28175(:,2));
    vlink_s6 = table2array(vlink_28175(:,3));
    vlink_s7 = table2array(vlink_28175(:,4));
    vlink_s8 = table2array(vlink_28175(:,5));

for i = 1:8
    if i < 5

        vlink_8145{:,i+1} = vlink_8145{:,i+1}-mean(vlink_8145{:,i+1},'omitnan');
    else
        vlink_28175{:,i-3} = vlink_28175{:,i-3}-mean(vlink_28175{:,i-3},'omitnan');
    end
end
%%
vel = movmean(ADV2_xvel,64,'omitnan');
figure
good = ADV2_xvel < -1 & ADV2_xvel > -3;
ax1 = subplot(211);
plot(ADV2_datetime_UTC,vel)
xlim([datetime('2022-11-22 19:49:00') datetime('2022-11-22 19:50:00')])
ylabel('Velocity [m/s]')
title('Blade Strain and Flow Velocity')
yyaxis right

ax = gca;
ax.YColor = 'k';
plot(vlink_8145{:,1},vlink_8145{:,4},'k')
ylabel('Microstrain')
ylim([-60 40])
grid on


ax2 = subplot(212);

plot(LC_datetime_UTC,Thrust_Force)
ylabel('Thrust [kN]')
yyaxis right
hold on
grid on
title('Blade Strain and Turbine Thrust')
 ylim([-60 40])
plot(vlink_8145{:,1},vlink_8145{:,4},'k')
ax = gca;
ax.YAxis(2).Color = 'k';
ylabel('Microstrain')
set(gcf,'color','w');
xlim([datetime('2022-11-22 19:49:00') datetime('2022-11-22 19:50:00')])
linkaxes([ax1 ax2],'x')

%% xy stuff

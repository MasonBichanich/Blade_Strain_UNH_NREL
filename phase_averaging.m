% load in strain data
A = exist('vlink_s1');

if A == 0
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
end

%% examine the data
plot(vlink_s1_datenum,vlink_s3)
title('Raw Data at S3')
zero = mean(vlink_s3);
hold on
yline(zero,'k')

%% chop up phases at peaks instead of at zero. The zero appears to have fluctuations over time, so I would need to do moving mean or something
% I've decided that the troughs of S3 are the most defined so I'll be using
% those
[dummy, peak_ind] = findpeaks(-vlink_s3,MinPeakDistance=100,MinPeakProminence=10);         %find the troughs as they are best defined
plot(vlink_s1_datenum(peak_ind),vlink_s3(peak_ind),'*r')                                     %confirm the peaks
title('Raw data truncated to turbine operation and highlighted peaks')

%% I'm going to see what happens if i remove the means from the data first
for i = 1:8
    if i < 5

        vlink_8145{:,i+1} = vlink_8145{:,i+1}-mean(vlink_8145{:,i+1},'omitnan');
    else
        vlink_28175{:,i-3} = vlink_28175{:,i-3}-mean(vlink_28175{:,i-3},'omitnan');
    end
end
%% break up the data at the troughs

%preallocation
broken = {cell(length(peak_ind),1), cell(length(peak_ind),1), cell(length(peak_ind),1), cell(length(peak_ind),1),...
    cell(length(peak_ind),1), cell(length(peak_ind),1), cell(length(peak_ind),1), cell(length(peak_ind),1),};
broken_time = cell(length(peak_ind),1);

%initializing
for i = 1:8
    if i < 5
        broken{i}{1} = vlink_28175(1:peak_ind(1),i+1);
    else
        broken{i}{1} = vlink_8145(1:peak_ind(1),i-3);
    end
end
broken_time{1} = vlink_s1_datenum(1:peak_ind(1));

%breaking the time
for i = 2:length(peak_ind)
       broken_time{i} = vlink_s1_datenum(peak_ind(i-1):peak_ind(i));
end

%breaking the data
for j = 1:8
    for i = 2:length(peak_ind)
        if j < 5
            broken{j}{i} = vlink_28175(peak_ind(i-1):peak_ind(i),j+1);
        else
            broken{j}{i} = vlink_8145(peak_ind(i-1):peak_ind(i),j-3);
        end
    end
end


%% check if it worked
for j = 1:8
    figure
    for i = 1:length(broken{1})
        plot(broken_time{i},table2array(broken{j}{i}))
        hold on
        title(num2str(j))
    end
end
%% figure for umerc
figure
for i = 1:length(broken{1})
    hold on
    plot(datetime(broken_time{i},'ConvertFrom','datenum'), ...
        table2array(broken{7}{i}))
end
ylabel('Microstrain')
set(gcf,'color','w');
xlim([datetime('2022-11-22 19:49:00') datetime('2022-11-22 19:50:00')])
xlabel('HH:MM:SS')
%% now lets average these mf phases. They are all different lengths. First try will be to make them all the same length 
% and Nan the missing data
% create a new compass, 0-360 and the length of each cell
% and use outerjoin

%create imaginary compass
for i = 1:height(broken{1})
    comp{i,:} = array2table(round((linspace(0,359,length(table2array(broken{1}{i})))))');
end

%combine imaginary compass with broken data
broken_comp = broken;
for j = 1:8
    for i = 1:height(broken{1})
        broken_comp{j}{i} = [comp{i},broken{j}{i}];
    end
end


%combining them all into a table (each degree has a bunch of nans and data
%points to average over) (this takes a long time)
tic
for j = 1:8
    joined{j} = (broken_comp{1}{1});
    for i = 2:height(broken_comp{1})
        joined{j} = outerjoin((broken_comp{j}{i}),joined{j},Keys=1,MergeKeys=true);
    end
    joined{j} = table2array(joined{j});
end
toc

%% the data should be ensemble averaged BEFORE being phase averaged.
dg = 5;
degree = num2str(dg);
ens_ave=[];
for k = 1:8
    for i = 1:length(joined{1})
        for j = 1:length(joined{1}(:,1))/dg
            ens_ave{k}{i}(j,:) = mean(joined{k}(dg*(j-1)+1:dg*j,i),'omitnan');
        end
    end
end

for i = 1:8
    ens_ave{i} = cell2mat(ens_ave{i});
end
%%
figure
cmap = jet(8)/1.2;
for i = 1:8

    phase_averaged{i} = [ens_ave{i}(:,1),mean(ens_ave{i}(:,2:end),2,'omitnan')];

    %subtract the mean
    means{i} = mean(phase_averaged{i}(:,2));
    pa_a{i} = phase_averaged{i}(:,2)-means{i};

    % shifting the zero crossing since its arbitrary
    pa_end = pa_a{i}(round(13/16*length(phase_averaged{i})):end); 
    pa_start = pa_a{i}(1:round(13/16*length(phase_averaged{i})-1));
    pa_shifted{i} = [pa_end;pa_start];
    
    hold on
    plot(phase_averaged{i}(:,1),phase_averaged{i}(:,2),'LineWidth',1,'Color', cmap(i,:))
end
yline(0)
legend('SG-1','SG-2','SG-3','SG-4','SG-5','SG-6','SG-7','SG-8')
grid on
xlim([0 360])
xlabel('Degree of rotation')
ylabel('Microstrain')
set(gcf,'color','w');
%%
figure
plot(phase_averaged{3}(:,1),phase_averaged{3}(:,2),'LineWidth',1,'Color', cmap(3,:))
yline(0)
legend('SG-3')
grid on
xlim([0 360])
xlabel('Degree of rotation')
ylabel('Microstrain')
set(gcf,'color','w');
%title("Phase average at sensor " + i + " with " + degree + " degree ensemble averages")

% figure
% plot(phase_averaged(:,1),phase_averaged(:,2))
% figure
% plot(linspace(-pi,pi,length(phase_averaged(:,1))),phase_averaged(:,2))

% pa_end = phase_averaged(round(1/4*length(phase_averaged)):end,2);
% pa_start = phase_averaged(1:round(1/4*length(phase_averaged)-1),2);
% pa_fixed = [pa_end;pa_start];
% figure
% plot(ens_ave(:,1),pa_fixed)
% grid on
% degree = num2str(dg);
% title("Phase averaged data with " + degree + " degree ensemble averages")

% %% averaging only the similar length phases
% j=1;
% for i = 1:1966
%     new_bc{i} = height(broken_comp{3}{i}) > 175 & height(broken_comp{3}{i}) < 300;
% end
% 
% tic
%     njoined = (broken_comp{3}{1});
%     for i = 2:height(broken_comp{1})
%         if new_bc{i} == 1
%             njoined = outerjoin((broken_comp{3}{i}),njoined,Keys=1,MergeKeys=true);
%         end
%     end
%     njoined = table2array(njoined);
% toc
% %%
% dg = 5;
% degree = num2str(dg);
% nens_ave = [];
% for i = 1:183
%     for j = 1:length(njoined(:,1))/dg
%         nens_ave{i}(j,:) = mean(njoined(dg*(j-1)+1:dg*j,i),'omitnan');
%     end
% end
% nens_ave = cell2mat(nens_ave);
% %%
%     nphase_averaged = [nens_ave(:,1),mean(nens_ave(:,2:end),2,'omitnan')];
% 
%     %subtract the mean
%     nmeans = mean(nphase_averaged(:,2),'omitnan');
%     npa_a = nphase_averaged(:,2)-nmeans;
% 
%     % shifting the zero crossing since its arbitrary
%     npa_end = npa_a(round(1/4.5*length(nphase_averaged)):end); 
%     npa_start = npa_a(1:round(1/4.5*length(nphase_averaged)-1));
%     npa_shifted = [npa_end;npa_start];
% 
%     hold on
%     plot(nphase_averaged(:,1),npa_shifted,'LineWidth',1)
% 
%% histogram of lengths of phases
Voltsys_table= load(['../../../Data_Files/MODAQ/Voltsys/Processed_Data/Daily_Files_No_QC/2022.11.22_Voltsys_No_QC_UTC_NREL_Final.mat']);
Voltsys_time_vector_UTC_Final_NREL = table2array(Voltsys_table.data(:,1:6));
round_sec =  0; % 1 if yes round to nearest second | 0 for no rounding
reverse   =  1; % 1 if you want to convert back to one cell timestamp  0 if you want to go to 6 cell timestamp (NREL Format)

[Voltsys_datetime_UTC] = f_newtime_NREL(Voltsys_time_vector_UTC_Final_NREL,round_sec,reverse);

Turbine_Freq = table2array(Voltsys_table.data(:,8));
Turbine_RPM = Turbine_Freq.*(120/40); % Speed of Turbine in RPM

Turbine_Spin_ind = (Turbine_RPM>1) & Voltsys_datetime_UTC > '2022-11-22 18:59:58.812' & Voltsys_datetime_UTC < '2022-11-22 20:27:24.546';

Turbine_RPM_avg = mean(Turbine_RPM(Turbine_Spin_ind));
RPM = Turbine_RPM(Turbine_Spin_ind);

Seconds_per_rotation_avg = (Turbine_RPM_avg/60)^-1;

dumpload_pwr_1s = table2array(Voltsys_table.data(:,12));

Capacitor_Voltage = table2array(Voltsys_table.data(:,48));

subplot(121)
histogram(RPM,60)

xlabel('Turbine RPM as measured by Voltsys')
ylabel('Counts')
annotation('textbox', [.32 .7 .1 .1],'String',"Mean = "+ mean(RPM) + newline + "Variance = "+var(RPM))
set(gcf,'color','w');
xlim([15 30])
for i = 1:length(broken_time{1})
    lens(i) = length(broken_time{i});
end
subplot(122)
histogram(64*60./lens,20)
ylabel('Counts')
xlabel('Turbine RPM as estimated by phase averaging algorithm')
annotation('textbox', [.77 .7 .1 .1],'String',"Mean = "+ mean(64*60./lens) + newline + "Variance = "+var(64*60./lens))
xlim([15 30])
set(gcf,'color','w');

Ztest = ((mean(RPM)-mean(64*60./lens))/sqrt(std(RPM)/length(RPM)+std(64*60./lens)/length(lens)));

[h,p] = ttest2(RPM,64*60./lens);

%%
%load('../../Data_Files/MODAQ/Load_Cells/Processed_Data/Daily_Files_Post_QC_no_timefix/2022.11.22_LC_Post_QC_UTC_NREL_Final.mat');
Load_Cell_table = load('../../../Data_Files/MODAQ/Load_Cells/Processed_Data/Daily_Files_Post_QC_no_timefix/2022.11.22_to_2022.11.22_LC_Post_QC_UTC_NREL_Final.mat');
% Load in Calibration Results

    Load_Cell_Cal_Results = load('../../../Data_Files/MODAQ/Load_Cells/Load_Cell_Calibration/Load_Cell_Cal_Results.mat');

% Load in Static Weight
%   Determined from LOAD_CELL_3_FORM_DRAG_ESTIMATE.m script

    load('../../../Data_Files/MODAQ/Load_Cells/Processed_Data/Turbine_Static_Weight_on_LC.mat')
     %% Extract Load Cell Variables    
LC_time_vector_UTC_Final_NREL = table2array(Load_Cell_table.data(:,1:6));
round_sec =  0; % 1 if yes round to nearest second | 0 for no rounding
reverse   =  1; % 1 if you want to convert back to one cell timestamp  0 if you want to go to 6 cell timestamp (NREL Format)

[LC_datetime_UTC] = f_newtime_NREL(LC_time_vector_UTC_Final_NREL,round_sec,reverse);

Port_LC = table2array(Load_Cell_table.data(:,7));
Star_LC = table2array(Load_Cell_table.data(:,8));

LC_Sum = Port_LC + Star_LC;
    b = New_static_weight; 
    
    m = -1/2*(abs(table2array(Load_Cell_Cal_Results.slopes(1,1))) + abs(table2array(Load_Cell_Cal_Results.slopes(1,3))));

    Thrust_Force = (LC_Sum - b)./m;

%%

plot(LC_datetime_UTC,Thrust_Force)
ylabel('Thrust [kN]')
yyaxis right
hold on
 
plot(vlink_8145{:,1},vlink_8145{:,4},'k')
ax = gca;
ax.YAxis(2).Color = 'k';
ylabel('Microstrain')
set(gcf,'color','w');
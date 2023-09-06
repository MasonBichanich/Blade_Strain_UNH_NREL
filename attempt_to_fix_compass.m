% The purpose of this script is to view the blade strain data on 11/22/22
% Output pertinent variables all on same clock?
close all
addpath('../../Functions')

A = exist("Voltsys_table");
if A == 0 
%% Load in the Data

% Load in the Load Cell data

%load('../../Data_Files/MODAQ/Load_Cells/Processed_Data/Daily_Files_Post_QC_no_timefix/2022.11.22_LC_Post_QC_UTC_NREL_Final.mat');
Load_Cell_table = load('../../Data_Files/MODAQ/Load_Cells/Processed_Data/Daily_Files_Post_QC_no_timefix/2022.11.22_to_2022.11.22_LC_Post_QC_UTC_NREL_Final.mat');
% Load in Calibration Results

    Load_Cell_Cal_Results = load('../../Data_Files/MODAQ/Load_Cells/Load_Cell_Calibration/Load_Cell_Cal_Results.mat');

% Load in Static Weight
%   Determined from LOAD_CELL_3_FORM_DRAG_ESTIMATE.m script

    load('../../Data_Files/MODAQ/Load_Cells/Processed_Data/Turbine_Static_Weight_on_LC.mat')
%%
% Load in Bow Side of Moon pool (ADV2) NO QC

ADV2_table = load('../../Data_Files/MODAQ/Vector_Velocity_Data_2/Processed_Data/2022_11_22_to_11_22_Vector_2_UTC.mat');
%%
% Load in the Voltsys Data

 Voltsys_table= load(['../../Data_Files/MODAQ/Voltsys/Processed_Data/Daily_Files_No_QC/2022.11.22_Voltsys_No_QC_UTC_NREL_Final.mat']);

 %% Load in the Vlink Data


 vlink_8145 = readtable('../../Data_Files/MODAQ/Vlink/Processed_Data/8145_vLink_strain-2022-11-22.csv');

 vlink_28175 = readtable('../../Data_Files/MODAQ/Vlink/Processed_Data/28175_vLink_strain-2022-11-22.csv');

 vlink_compass = readtable('../../Data_Files/MODAQ/Vlink/Processed_Data/compass_vLink-2022-11-22.csv');


  
 %% Vlink Positional Data 
 %  (listed from blade top nearest free surface when deployed) to bottom)
% FB01-10CM     VLink 1 - SN 8145	1	1
% FB02-23CM     VLink 1 - SN 8145	1	2
% FB03-36CM     VLink 1 - SN 8145	1	3
% FB04-70CM     VLink 1 - SN 8145	1	4
% FB05-85CM     VLink 2 - SN 28175	2	1
% FB06-100CM	VLink 2 - SN 28175	2	2
% FB07-134CM	VLink 2 - SN 28175	2	3
% FB08-160CM	VLink 2 - SN 28175	2	4
 %% Extract Compass Variables    
  
vlink_compass_time = table2array(vlink_compass(:,1));
vlink_compass_data = table2array(vlink_compass(:,2));

% Extract vlink 8145 data
vlink_s1_s4_time = table2array(vlink_8145(:,1));
vlink_s1 = table2array(vlink_8145(:,2));
vlink_s2 = table2array(vlink_8145(:,3));
vlink_s3 = table2array(vlink_8145(:,4));
vlink_s4 = table2array(vlink_8145(:,5));

% Extract vlink 28175 data
vlink_s5_s8_time = table2array(vlink_28175(:,1));
vlink_s5 = table2array(vlink_28175(:,2));
vlink_s6 = table2array(vlink_28175(:,3));
vlink_s7 = table2array(vlink_28175(:,4));
vlink_s8 = table2array(vlink_28175(:,5));

 %% Extract Load Cell Variables    
LC_time_vector_UTC_Final_NREL = table2array(Load_Cell_table.data(:,1:6));
round_sec =  0; % 1 if yes round to nearest second | 0 for no rounding
reverse   =  1; % 1 if you want to convert back to one cell timestamp  0 if you want to go to 6 cell timestamp (NREL Format)

[LC_datetime_UTC] = f_newtime_NREL(LC_time_vector_UTC_Final_NREL,round_sec,reverse);

Port_LC = table2array(Load_Cell_table.data(:,7));
Star_LC = table2array(Load_Cell_table.data(:,8));

LC_Sum = Port_LC + Star_LC;

%% Extract ADV2 Variable

ADV2_datenum_UTC  = table2array(ADV2_table.Vector_new(:,1));

ADV2_datetime_UTC = datetime(ADV2_datenum_UTC,'ConvertFrom','datenum','Format', 'MM/dd/yyyy HH:mm:ss.SSSSSSSSS');

ADV2_xvel = table2array(ADV2_table.Vector_new(:,2));

figure
    plot(ADV2_datetime_UTC,ADV2_xvel,'b .')

%% Extract Voltsys Variables    
Voltsys_time_vector_UTC_Final_NREL = table2array(Voltsys_table.data(:,1:6));
round_sec =  0; % 1 if yes round to nearest second | 0 for no rounding
reverse   =  1; % 1 if you want to convert back to one cell timestamp  0 if you want to go to 6 cell timestamp (NREL Format)

[Voltsys_datetime_UTC] = f_newtime_NREL(Voltsys_time_vector_UTC_Final_NREL,round_sec,reverse);

Turbine_Freq = table2array(Voltsys_table.data(:,8));
Turbine_RPM = Turbine_Freq.*(120/40); % Speed of Turbine in RPM

Turbine_Spin_Indices = find(Turbine_RPM~=0);
Turbine_RPM_avg = mean(Turbine_RPM(Turbine_Spin_Indices));

Seconds_per_rotation_avg = (Turbine_RPM_avg/60)^-1;

dumpload_pwr_1s = table2array(Voltsys_table.data(:,12));

Capacitor_Voltage = table2array(Voltsys_table.data(:,48));

%% Convert to Turbine Thrust Force

% Convert LC results to Drag Force Estimates using LoadCell Calibration Results 
% (Load_Cell_Cal_Results structure)

    % NEED TO CHECK ON if this is the proper method of selecting the y intercept
    b = New_static_weight; 

    m = -1/2*(abs(table2array(Load_Cell_Cal_Results.slopes(1,1))) + abs(table2array(Load_Cell_Cal_Results.slopes(1,3))));

    Thrust_Force = (LC_Sum - b)./m;
else
end

%% Convert ADV data to horizontal velocity magnitude
% Run ADV data through QC script
% Update Voltsys output data to reflect changes in dumpload power levels
% (./10)

%% get rid of that weird resampling thing thats going on. Ask about this
jumps = diff(vlink_compass_data);
not_flat_ind = abs(jumps) > 5;
figure
plot(vlink_compass_time(not_flat_ind),vlink_compass_data(not_flat_ind))
new_comp = vlink_compass_data(not_flat_ind);
new_comp_time = vlink_compass_time(not_flat_ind);
hold on
plot(vlink_compass_time,vlink_compass_data)
%plot(vlink_compass_time(2:end),compass_data_delta,'b *')
   date_begin = datetime('2022,11,22,19,55,00','InputFormat','yyyy,MM,dd,HH,mm,ss','format','MM/dd/yyyy HH:mm:ss');
   date_end = datetime('2022,11,22,20,00,30','InputFormat','yyyy,MM,dd,HH,mm,ss','format','MM/dd/yyyy HH:mm:ss');
   xlim([date_begin date_end])
%%     
    %compass_data_jump_indices = find(compass_data_delta<-200);    %What should this number be? a better way of doing this may be local minima
    [dummy, compass_data_jump_indices] = findpeaks(new_comp,'MinPeakProminence',200);

    vlink_compass_time_jump = new_comp_time(compass_data_jump_indices);
    figure
    hold on
    plot(new_comp_time,new_comp)
    plot(vlink_compass_time_jump,new_comp(compass_data_jump_indices),'r .')
    xlim([date_begin date_end])
%% Trial run for correcting compass 
indicies_per_rev =  diff(compass_data_jump_indices);
number_of_revs = length(compass_data_jump_indices);
n = 0;  
for i = 1:length(new_comp)    
    seq_comp(i) = new_comp(i)+(n)*360;                          %This loop makes the compass sequential instead of reseting at 360
    if any(compass_data_jump_indices == i)
        n = n + 1;
    end
end
% seq_comp = seq_comp';
% seq_comp_jump = abs(diff(seq_comp)) > 20;
% seq_comp(seq_comp_jump) = NaN;

%seq_comp_corrected = seq_comp - 2;          %number of degrees of compass drift

B = exist('ratio');
if B == 1
    seq_comp_corrected = seq_comp*ratio;          %scale factor for compass based on total RPM measured by voltsys
else
    seq_comp_corrected = seq_comp;
end

figure

plot(new_comp_time,seq_comp)
xlim([date_begin date_end])
%Now we need to go back to regular compass 
n = 0;
for i = 1:length(new_comp)
    if seq_comp_corrected(i) > (n+1)*360              
       n = n + 1;
    end
    corrected_comp(i) = seq_comp_corrected(i)-n*360;

end

corrected_comp(1) = 0;
corrected_comp = corrected_comp';
figure
plot(new_comp_time,corrected_comp)
ylim([0 360])
xlim([date_begin date_end])
%%     have to repeat this section with the new time series as the peaks have changed
    %compass_data_jump_indices = find(compass_data_delta<-200);    %What should this number be? a better way of doing this may be local minima
    [dummy, compass_data_jump_indices_corr] = findpeaks(corrected_comp,'MinPeakProminence',200);

    vlink_compass_time_jump = new_comp_time(compass_data_jump_indices_corr);
    figure
    hold on
    plot(new_comp_time,corrected_comp)
    plot(vlink_compass_time_jump,corrected_comp(compass_data_jump_indices_corr),'r .')

    %% Now determine the time of each rotation
    
    % for i = 1:length(vlink_compass_time_jump)-1
    %     compass_1rev_times(i) = vlink_compass_time_jump(i+1)-vlink_compass_time_jump(i);
    % end
    compass_1rev_times = diff(vlink_compass_time_jump);
    
    figure
    plot(compass_1rev_times,' b.--')
    seconds_compass_1rev_times = seconds(compass_1rev_times);
    
    
%% Select data indices that occurs before each good 1 rev time < 3.5 s

    good_compass_1rev_times = find(seconds_compass_1rev_times<4 & seconds_compass_1rev_times>1.7);

    good_rev_times = vlink_compass_time_jump(good_compass_1rev_times);
    
    % the original vector indices that coincide with these good times
     good_rev_indices = compass_data_jump_indices_corr(good_compass_1rev_times);
    
    % Take the time of each good rev and determine the number of data
    % points that should exist per rev based on sample rate
    delta_t_compass = 1/10;
    good_rev_indices_min = round(good_rev_indices - seconds_compass_1rev_times(good_compass_1rev_times)./delta_t_compass,0);
    
    % figure
    %     plot(new_comp_time,new_comp,'k --')
    % 
        % for i =2:length(good_rev_indices_min)
        % hold on
        % plot(new_comp_time(good_rev_indices_min(i):good_rev_indices(i)),new_comp(good_rev_indices_min(i):good_rev_indices(i)),'r*--')
        % end
        
    %% How many indices are in each good rev range
    
    good_rev_indices_length = good_rev_indices - good_rev_indices_min;
    
    figure
    plot(good_rev_indices_length,'k *')
    
    max(good_rev_indices_length);
    
    min_index =  find(min(good_rev_indices_length));
    
    %% need to scale each rev by the number of points
    good_rev_indices_length_scale = good_rev_indices_length./min(good_rev_indices_length);

   %good_rev_indices_length_scale = good_rev_indices_length_scale(:,1);

    %% try selecting indices within a specific range and adjusted by scaling for each rev
    clear trial_output
    
    % Which portion of each revolution do you want to view
    prompt = 'Which portion of each revolution do you want to view lower limit[as a % of 1 revolution]? ';
    Screen1 = input(prompt);
    
    prompt = 'Which portion of each revolution do you want to view upper limit[as a % of 1 revolution]? ';
    Screen2 = input(prompt);
        
    
     trial_range = [Screen1/100 Screen2/100]; % ie indices ranging from 20% to 30% over each rpm
    trial_range_string = strcat(num2str(trial_range(1)*100),"% to ",num2str(trial_range(2)*100),"% of one revolution");
    
    % Which sensor do you want to view
    prompt = 'What Sensor do you want to plot [s1 =1 s2 = 2 ... s8 == 8]? ';
    Screen = input(prompt);

    if (Screen == 1)
            k =1;
            sensor_ID = "S1";
      elseif(Screen == 2)
            k =2;
            sensor_ID = "S2";
      elseif(Screen == 3)
            k =3;
            sensor_ID = "S3";
      elseif(Screen == 4)
            k =4;
            sensor_ID = "S4";
      elseif(Screen == 5)
            k =5;
            sensor_ID = "S5";
      elseif(Screen == 6)
            k =6;        
            sensor_ID = "S6";
      elseif(Screen == 7)
            k =7;
            sensor_ID = "S7";
      elseif(Screen == 8)
            k =8;
            sensor_ID = "S8";
    end    
     figure
    %ax1 = subplot(211);
        if k == 1
            plot(vlink_s1_s4_time,vlink_s1,'k .--')
        elseif k ==2
            plot(vlink_s1_s4_time,vlink_s2,'k .--')
        elseif k ==3
            plot(vlink_s1_s4_time,vlink_s3,'k .--')    
        elseif k ==4
            plot(vlink_s1_s4_time,vlink_s4,'k .--')    
        elseif k ==5
            plot(vlink_s5_s8_time,vlink_s5,'k .--')    
        elseif k ==6
            plot(vlink_s5_s8_time,vlink_s6,'k .--') 
        elseif k ==7
            plot(vlink_s5_s8_time,vlink_s7,'k .--')     
        elseif k ==8
            plot(vlink_s5_s8_time,vlink_s8,'k .--') 
        end
    hold on
    %yyaxis right
    % plot(vlink_compass_time,vlink_compass_data,' b --')   
    
%       hold on
%        legend('s1','compass','s1 range','compass range')

    % compass measurements within range
    hold on
    %yyaxis right
    comp_ind = corrected_comp > 360*Screen1/100 & corrected_comp < 360*Screen2/100;
    %plot(new_comp_time(comp_ind),corrected_comp(comp_ind),'*g')
    time_in_range = new_comp_time(comp_ind);
    time_ranges = (ismember(vlink_s1_s4_time,time_in_range));
    % for i = 1:length(vlink_s1_s4_time)
    %     time_range(i) = any(time_in_range == vlink_s1_s4_time(i));
    % end
    %yyaxis left
        if k == 1
            plot(vlink_s1_s4_time(time_ranges),vlink_s1(time_ranges),'*r')
        elseif k ==2
            plot(vlink_s1_s4_time(time_ranges),vlink_s2(time_ranges),'*r')
        elseif k ==3
            plot(vlink_s1_s4_time(time_ranges),vlink_s3(time_ranges),'*r')
        elseif k ==4
            plot(vlink_s1_s4_time(time_ranges),vlink_s4(time_ranges),'*r')
        elseif k ==5
            plot(vlink_s5_s8_time(time_ranges),vlink_s5(time_ranges),'*r')
        elseif k ==6
            plot(vlink_s5_s8_time(time_ranges),vlink_s6(time_ranges),'*r')
        elseif k ==7
            plot(vlink_s5_s8_time(time_ranges),vlink_s7(time_ranges),'*r')
        elseif k ==8
            plot(vlink_s5_s8_time(time_ranges),vlink_s8(time_ranges),'*r')
        end
    
    % for trial_index = 2:length(good_rev_indices_length_scale)
    % 
    % trial_output(trial_index).range = [round(good_rev_indices_length_scale(trial_index)*trial_range(1)*min(good_rev_indices_length) ) ...
    %                         +      good_rev_indices_min(trial_index), ...
    %                       round(good_rev_indices_length_scale(trial_index)*trial_range(2)*min(good_rev_indices_length) ) ...
    %                         +      good_rev_indices_min(trial_index)
    %                       ];
    % yyaxis left
    % 
    %     if k == 1
    %     %plot(vlink_s1_s4_time(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),vlink_s1(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),'r *')
    %     %plot(vlink_s1_s4_time(time_ranges),vlink_s4(time_ranges),'*r')
    %     elseif k ==2
    %     plot(vlink_s1_s4_time(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),vlink_s2(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),'r *')
    %     elseif k ==3
    %     plot(vlink_s1_s4_time(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),vlink_s3(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),'r *')
    %     elseif k ==4
    %     plot(vlink_s1_s4_time(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),vlink_s4(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),'r *')
    %     elseif k ==5
    %     plot(vlink_s5_s8_time(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),vlink_s5(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),'r *')
    %     elseif k ==6
    %     plot(vlink_s5_s8_time(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),vlink_s6(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),'r *')
    %     elseif k ==7
    %     plot(vlink_s5_s8_time(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),vlink_s7(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),'r *')
    %     elseif k ==8
    %     plot(vlink_s5_s8_time(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),vlink_s8(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),'r *')
    %     end
    % 
    % % hold on
    % % yyaxis right
    % %plot(new_comp_time(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),corrected_comp(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),' g *')  
    % hold on
    % end  
    % 

  %ax2 = subplot(212);
  %plot(ADV2_datetime_UTC,ADV2_xvel.^3,'-.')
  %ylim([-9 2])
%linkaxes([ax1 ax2],'x')
   date_begin = datetime('2022,11,22,19,00,00','InputFormat','yyyy,MM,dd,HH,mm,ss','format','MM/dd/yyyy HH:mm:ss');
   date_end = datetime('2022,11,22,20,25,30','InputFormat','yyyy,MM,dd,HH,mm,ss','format','MM/dd/yyyy HH:mm:ss');
    
   xlim([date_begin date_end])
   
   title_string = [sensor_ID,trial_range_string];
   title(title_string)
%%
time = minutes(vlink_compass_time_jump(end)-vlink_compass_time_jump(1));
avg_RPM = mean(Turbine_RPM(Voltsys_datetime_UTC > vlink_compass_time_jump(1) & Voltsys_datetime_UTC < vlink_compass_time_jump(end)));
revs = avg_RPM*time;
voltsys_degrees = revs*360;
last_value = seq_comp_corrected(end);
% compass_revs = length(compass_data_jump_indices);
correction = last_value - voltsys_degrees;
% num_to_add = correction*360;

%new_last_value = last_value-num_to_add;
C = exist('ratio');
if C == 0
    ratio = voltsys_degrees/last_value
end
figure
yyaxis left; 
plot(vlink_s1_s4_time,vlink_s1,'k'); 
hold on; 
yyaxis right; 
plot(new_comp_time,corrected_comp,'.--b');

%% spectra of strain

[S_v,f] = f_spectra(vlink_s1,1/64,5,5);  %2187000:2240000
title('Spectrum of S1')
[S_c,f] = f_spectra(new_comp,1/10,5,5);
title('Spectrum of compass')

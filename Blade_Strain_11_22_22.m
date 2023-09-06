% The purpose of this script is to view the blade strain data on 11/22/22
% Output pertinent variables all on same clock?
addpath('../../Functions')

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
 current_path = pwd;
 
 cd('../../Data_Files/MODAQ/Vlink/Processed_Data')
 
 vlink_8145 = readtable('8145_vLink_strain-2022-11-22.csv');
 
 vlink_28175 = readtable('28175_vLink_strain-2022-11-22.csv');
 
 vlink_compass = readtable('compass_vLink-2022-11-22.csv');

 cd(current_path)
  
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

%% Convert ADV data to horizontal velocity magnitude
% Run ADV data through QC script
% Update Voltsys output data to reflect changes in dumpload power levels
% (./10)

%% All Strain gauges plot

figure (7)
    plot(vlink_s1_s4_time,vlink_s1,'c-')
    hold on
    plot(vlink_s1_s4_time,vlink_s2,'k -')
    hold on    
    plot(vlink_s1_s4_time,vlink_s3,'g-')
    hold on
    plot(vlink_s1_s4_time,vlink_s4,'r -')
    
    hold on
    
    plot(vlink_s5_s8_time,vlink_s5,'y --')
    hold on    
    plot(vlink_s5_s8_time,vlink_s6,'r--')
    hold on      
    plot(vlink_s5_s8_time,vlink_s7,'g --')
    hold on      
    plot(vlink_s5_s8_time,vlink_s8,'c --')
 
    yline (0)
    ylim([-100 100])
    legend('s1','s2','s3','s4', 's5','s6','s7','s8')
    xlabel('Time')
    ylabel('Strain')
figure (7)
hold on
    yyaxis right
    plot(vlink_compass_time,vlink_compass_data,'k.-')
    
%% select indicies corresponding to compass position range 
% signal roughly 0 to 360 in value
clear indices_90
indices_90 = find(vlink_compass_data>85 & vlink_compass_data<95);

figure(7)
hold on
    yyaxis right
    hold on
    plot(vlink_compass_time(indices_90),vlink_compass_data(indices_90),'r.-')

figure(7)
hold on
    yyaxis left
    plot(vlink_s1_s4_time(indices_90),vlink_s1(indices_90),'c *')
    hold on
    plot(vlink_s1_s4_time(indices_90),vlink_s2(indices_90),'k *')
    hold on    
    plot(vlink_s1_s4_time(indices_90),vlink_s3(indices_90),'g *')
    hold on
    plot(vlink_s1_s4_time(indices_90),vlink_s4(indices_90),'r *')
    
    hold on
    
    plot(vlink_s5_s8_time(indices_90),vlink_s5(indices_90),'y *')
    hold on    
    plot(vlink_s5_s8_time(indices_90),vlink_s6(indices_90),'r *')
    hold on      
    plot(vlink_s5_s8_time(indices_90),vlink_s7(indices_90),'g *')
    hold on      
    plot(vlink_s5_s8_time(indices_90),vlink_s8(indices_90),'c *')
    
%% One Strain Gage
figure(8)
    yyaxis right
    plot(vlink_compass_time,vlink_compass_data,'k-')
    hold on
    plot(vlink_compass_time(indices_90),vlink_compass_data(indices_90),'r*')

figure(8)
hold on
    yyaxis left
    plot(vlink_s1_s4_time,vlink_s1,'c-')
    hold on
    plot(vlink_s1_s4_time(indices_90),vlink_s1(indices_90),'c *')
    hold on

    %%
%  figure
%     scatter(vlink_compass_data(indices_90),vlink_s1(indices_90),'k *')
%     
%     hold on
% 
%     scatter(vlink_compass_data(indices_90(1:100)),vlink_s1(indices_90(1:100)),'r .')

    vlink_s1_s4_time(indices_90(1));
    
    vlink_s1_s4_time(indices_90(100));

    delta_t_compass = seconds(vlink_s1_s4_time(100)-vlink_s1_s4_time(99));
    
    compass_samples_per_rev = Seconds_per_rotation_avg/delta_t_compass;
    
    %% So roughly every 171 samples of position per revolution
    % Compare strain readings every 171 points from one another
    
    compare_indices = 1:171:length(vlink_s5_s8_time);
    
    figure
    plot(vlink_s1_s4_time,vlink_s1,'k .-')
    hold on
    plot(vlink_s1_s4_time(compare_indices),vlink_s1(compare_indices),'r *')
    hold on
    yyaxis right
    hold on
    plot(vlink_compass_time,vlink_compass_data,'b .-')
    plot(vlink_compass_time(compare_indices),vlink_compass_data(compare_indices),'g *')

    %% Find compass indices when it resets
    % by looking for large delta in compass value according with wrapping
    % from 360 back to 0
    
    % for i = 1:length(vlink_compass_data)-1
    % compass_data_delta(i) = vlink_compass_data(i+1)-vlink_compass_data(i);
    % end
    
    figure
    %plot(vlink_compass_time(2:end),compass_data_delta,'b *')
     
    %compass_data_jump_indices = find(compass_data_delta<-200);    %What should this number be? a better way of doing this may be local minima
    [dummy, compass_data_jump_indices] = findpeaks(vlink_compass_data,'MinPeakProminence',200,'MinPeakDistance',10);
    %vlink_compass_time_jump = vlink_compass_time(compass_data_jump_indices);
    vlink_compass_time_jump = vlink_compass_time(compass_data_jump_indices);
    
    hold on
    plot(vlink_compass_time,vlink_compass_data)
    plot(vlink_compass_time_jump,vlink_compass_data(compass_data_jump_indices),'r .')
    
    
    %% Now determine the time of each rotation
    
    % for i = 1:length(vlink_compass_time_jump)-1
    %     compass_1rev_times(i) = vlink_compass_time_jump(i+1)-vlink_compass_time_jump(i);
    % end
    compass_1rev_times = diff(vlink_compass_time_jump);
    
    figure
    plot(compass_1rev_times,' b.--')
    seconds_compass_1rev_times = seconds(compass_1rev_times);
    
    
%% Select data indices that occurs before each good 1 rev time < 3.5 s

    good_compass_1rev_times = find(seconds_compass_1rev_times<4);

    good_rev_times = vlink_compass_time_jump(good_compass_1rev_times);
    
    % the original vector indices that coincide with these good times
     good_rev_indices = compass_data_jump_indices(good_compass_1rev_times);
    
    % Take the time of each good rev and determine the number of data
    % points that should exist per rev based on sample rate
    
    good_rev_indices_min = round(good_rev_indices - seconds_compass_1rev_times(good_compass_1rev_times)./delta_t_compass,0);
    
    figure
        plot(vlink_compass_time,vlink_compass_data,'k --')
        
        for i =2:length(good_rev_indices_min)
        hold on
        plot(vlink_compass_time(good_rev_indices_min(i):good_rev_indices(i)),vlink_compass_data(good_rev_indices_min(i):good_rev_indices(i)),'r*--')
        end
        
    %% How many indices are in each good rev range
    
    good_rev_indices_length = good_rev_indices- good_rev_indices_min;
    
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
        if k == 1
            plot(vlink_s1_s4_time,vlink_s2,'k .--')
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
    yyaxis right
    plot(vlink_compass_time,vlink_compass_data,' b .--')   
    
%       hold on
%        legend('s1','compass','s1 range','compass range')

    for trial_index = 2:1000%length(good_rev_indices_length_scale)
    
    trial_output(trial_index).range = [round(good_rev_indices_length_scale(trial_index)*trial_range(1)*min(good_rev_indices_length) ) ...
                            +      good_rev_indices_min(trial_index) ...
                          round(good_rev_indices_length_scale(trial_index)*trial_range(2)*min(good_rev_indices_length) ) ...
                            +      good_rev_indices_min(trial_index)
                          ];
    yyaxis left
    
        if k == 1
        plot(vlink_s1_s4_time(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),vlink_s1(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),'r *')
        elseif k ==2
        plot(vlink_s1_s4_time(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),vlink_s2(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),'r *')
        elseif k ==3
        plot(vlink_s1_s4_time(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),vlink_s3(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),'r *')
        elseif k ==4
        plot(vlink_s1_s4_time(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),vlink_s4(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),'r *')
        elseif k ==5
        plot(vlink_s5_s8_time(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),vlink_s5(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),'r *')
        elseif k ==6
        plot(vlink_s5_s8_time(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),vlink_s6(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),'r *')
        elseif k ==7
        plot(vlink_s5_s8_time(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),vlink_s7(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),'r *')
        elseif k ==8
        plot(vlink_s5_s8_time(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),vlink_s8(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),'r *')
        end
           
    hold on
    yyaxis right
    plot(vlink_compass_time(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),vlink_compass_data(trial_output(trial_index).range(1):trial_output(trial_index).range(2)),' g *')   
    hold on
    end                      
    

   date_begin = datetime('2022,11,22,19,17,00','InputFormat','yyyy,MM,dd,HH,mm,ss','format','MM/dd/yyyy HH:mm:ss');
   date_end = datetime('2022,11,22,19,19,30','InputFormat','yyyy,MM,dd,HH,mm,ss','format','MM/dd/yyyy HH:mm:ss');
    
   xlim([date_begin date_end])
   
   title_string = [sensor_ID,trial_range_string];
   title(title_string)  
%% Test plot    
figure
tiledlayout(5,1)

ax1 = nexttile;
    plot(LC_datetime_UTC,Thrust_Force,'b.--')
    title('Thrust Load (kN)')
ax2 =nexttile
    plot(ADV2_datetime_UTC, ADV2_xvel,'b .--')
    title('ADV2 vel x (m/s)')
ax3 = nexttile;
    plot(Voltsys_datetime_UTC,dumpload_pwr_1s./10,'b .--')
    title('Dumpload 1s power (kW)')
   
ax4 = nexttile;
    plot(Voltsys_datetime_UTC,Turbine_RPM,'b .--')
    title('Turbine Shaft Speed (rpm)')

ax5 = nexttile;
     plot(table2array(vlink_compass(:,1)),table2array(vlink_compass(:,2)),'b .--')
     hold on
     yyaxis right
      plot(table2array(vlink_8145(:,1)),table2array(vlink_8145(:,4)),'r .--')  
      
linkaxes([ax1 ax2 ax3 ax4 ax5],'x')

%%
   figure
 plot(vlink_compass_time,vlink_compass_data,'b .--')
 hold on
 
 figure
yyaxis left
    plot(Voltsys_datetime_UTC,Turbine_RPM,'r .--')
 hold on
 yyaxis right
  plot(vlink_s1_s4_time,vlink_s1,'k .--')



    %%

%% Save Load Cell Table with Turbine Thrust Force Column


New_Load_Cell(:,(1:6)) = Load_Cell_table.data(:,(1:6)) ;
New_Load_Cell(:,(7:8)) = Load_Cell_table.data(:,(7:8)) ;
New_Load_Cell(:,9) = array2table(Thrust_Force(:,1)) ;

New_Load_Cell.Properties.VariableNames(1:8) = Load_Cell_table.data.Properties.VariableNames(1:8);
New_Load_Cell.Properties.VariableUnits(1:8) = Load_Cell_table.data.Properties.VariableUnits(1:8);


New_Load_Cell.Properties.VariableNames(9) = {'Turbine Thrust Force'};
New_Load_Cell.Properties.VariableUnits(9) = {'kN'};

filename = '2022.11.22_LC_Post_QC_UTC_NREL_Final_with_Thrust_Force.mat';

%save(filename,'New_Load_Cell','-v7.3')

% Save as a CSV File

filename = '2022.11.22_LC_Post_QC_UTC_NREL_Final_with_Thrust_Force.csv';

% write(New_Load_Cell,filename,'Delimiter',',')

units = New_Load_Cell.Properties.VariableUnits;
units = cell2table(units);

%write(units,'2022.11.22_to_2022.11.22_Load_Cell_Data__with_Thrust_Force_Variable_Units.csv','Delimiter',',')






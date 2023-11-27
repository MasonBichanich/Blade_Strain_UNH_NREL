%  The purpose of this script is to perform quality control tests on each
%  ADV independently

% Add Functions Driectory to Path
addpath('..\Blade_Strain_UNH_NREL\Functions')
%% Load each ADV seperately and then close all and repeat for the next instrument

%load in the data
clear
tic 
sf = 64;
%%   Which ADV data did you load in?
prompt = 'Which ADV? [1, 2, or 3]';
 str = input(prompt,'s');
tic

if (str2double(str) == 1)
    m = ' ADV 1 ';
    k =1;
elseif (str2double(str) == 2)
    m = ' ADV 2 ';
    k =2;
else 
    m = ' ADV 3 ';
    k = 3;
end   

% Initialize variable for ADV1 or ADV2 specific criteria to compare against k
      k1 = 1;
      k2 = 2;
      k3 = 3;
% %%  
% % BOW ADV 2 Data corrected to NOAA Tidal Direction Standard
% load("../../ADV-Nortek Vector/8.17.23 Deployment/ADV"+str+"/Converted Data/Vector_"+str+"_UTC_mod.mat")
%     Vector = Vector_new;
%     Vector_datetime_UTC = datetime(table2array(Vector(:,1)),'ConvertFrom','datenum','Format', 'MM/dd/yyyy HH:mm:ss.SSSSSSSSS');
% 
%     Vector_datenum_UTC = table2array(Vector(:,1));
% 
%     clear Vector_new
%     Vector_nums = table2array(Vector(:,2:13));
% 
% % BOW ADV 2 System Data
% % 
% load("../../ADV-Nortek Vector/8.17.23 Deployment/ADV"+str+"/Converted Data/Vector_"+str+"_SystemData_UTC.mat")
%     Vector_System_Data = Vector_SystemData_new;
%     Vector_System_Data_datetime_UTC = datetime(table2array(Vector_System_Data(:,1)),'ConvertFrom','datenum','Format', 'MM/dd/yyyy HH:mm:ss.SSSSSSSSS');
%     Vector_System_Data_datenum_UTC = table2array(Vector_System_Data(:,1));
%     clear Vector_SystemData_new
% 
% 
% 
% toc
%%
load('../../../Data_Files/MODAQ/Vector_Velocity_Data_2/Processed_Data/2022_11_22_to_11_22_Vector_2_UTC.mat')
    Vector = Vector_new;
    Vector_datetime_UTC = datetime(table2array(Vector(:,1)),'ConvertFrom','datenum','Format', 'MM/dd/yyyy HH:mm:ss.SSSSSSSSS');

    Vector_datenum_UTC = table2array(Vector(:,1));

    clear Vector_new
    Vector_nums = table2array(Vector(:,2:13));


%% Plot?
prompt = 'Would you like to display figures (Yes ==1/N0 ==0)? ';
PLOT = input(prompt);

if (prompt == 1)
    PLOT = input(prompt);
elseif (prompt == 0)
    PLOT = input(prompt);
end
%% Inititialize the test criteria and flag matrizes for each test
%
% Nortek Comprehensive Manual Page 38 for data output from Vectors
  
% Parameters to QC % Battery Voltage? % Sound Speed? % Temperature?
 
%Test 1
%     Battery_min = 21.5;  % Minimum acceptible battery voltage level
%     Battery_max = NaN;
%     
%     ADCP2_battery_nums = table2array(ADCP2.ADCP2(:,3));
% 
%     Flagged_Struct.btvlt = ones(height(ADCP2_battery_nums),width(ADCP2_battery_nums));
    
% Test 6 - Signal Strength Amplitude Beam 1 Beam 2 Beam 3
    % The scaling for the conversion from counts to dB varies a little from one instrument to another, 
    % but 1 count is around 0.4-0.45 dB. For more details check out the technical note (No. 003) 
    % named Monitoring Sediment Concentration with acoustic backscattering instruments available at the 
    %Nortek site (http://www.nortek-as.com/en/knowledge-center/technical-notes).

    Amp_2_db = 0.4;
    
if k == k2 
      SIGSTR_max = 250*Amp_2_db; % = 100 db
      SIGSTR_min = 45*Amp_2_db;  % = 18 db
else 
      SIGSTR_max = 250*Amp_2_db; % = 100 db
      SIGSTR_min = 51*Amp_2_db;  % = 20.4 db
end
      
%   SIGSTR_max = 160*Amp_2_db;
%   SIGSTR_min = 80*Amp_2_db;

    ADV_Amp1_nums = Vector_nums(:,6).*Amp_2_db;
    ADV_Amp2_nums = Vector_nums(:,7).*Amp_2_db;
    ADV_Amp3_nums = Vector_nums(:,8).*Amp_2_db;

    Flagged_Struct.SIGSTR = ones(height(ADV_Amp1_nums),width(ADV_Amp1_nums));

% Test 7 - Correlation Beam 1 Beam 2 Beam 3
    Corr_min = 70;
    %Corr_min = 0.3+ 0.4*sqrt(sf/25);
    Corr_max = 100;
        
    ADV_Corr1_nums = Vector_nums(:,9) ;
    ADV_Corr2_nums = Vector_nums(:,10) ;
    ADV_Corr3_nums = Vector_nums(:,11) ;

    Flagged_Struct.Corr = ones(height(ADV_Corr1_nums),width(ADV_Corr1_nums));
              
 % Test 10 - Current Magnitude
    ADV_Current_Magnitude_min = 0;
    ADV_Current_Magnitude_max = 5.25;

    Flagged_Struct.Horz_Current_Mag = ones( height(ADV_Corr1_nums), width(ADV_Corr1_nums));
           
% Test 11 - Current Direction
    ADV_Current_Direction_min = -180;
    ADV_Current_Direction_max = 180;

    Flagged_Struct.Horz_Current_Dir = ones( height(ADV_Corr1_nums), width(ADV_Corr1_nums));
            
% Test 12a - X Velocity Max Value
    % Horizontal Velocity - Instrument Coordinate Velocity X direction
    INSVELX_max_min = -4;
    INSVELX_max_max = 4;

    INSVELX_max_nums = Vector_nums(:,1);  

    Flagged_Struct.INSVELX_max = ones(height(INSVELX_max_nums),width(INSVELX_max_nums));

% Test 12b - Y Velocity Max Value
    % Horizontal Velocity - Instrument Coordinate Velocity
    INSVELY_max_min = -4;
    INSVELY_max_max = 4;

    INSVELY_max_nums = Vector_nums(:,2);  

    Flagged_Struct.INSVELY_max = ones(height(INSVELY_max_nums),width(INSVELY_max_nums));

% Test 13 - Z Velocity Max Value
    % Vertical Velocity - Instrument Coordinate Velocity 
    INSVELZ_max_min = -1.5;
    INSVELZ_max_max = 1.5;

    INSVELZ_max_nums = Vector_nums(:,3);  

    Flagged_Struct.INSVELZ_max = ones(height(INSVELZ_max_nums),width(INSVELZ_max_nums));

    
% Test 15a - Rate of Change XVEL
    INSVELX_dt_max = 1;
    INSVELX_dt_min = NaN;
    
    ADV_INSVELX_dt_nums =  Vector_nums(:,1);

    Flagged_Struct.INSVELX_dt = ones(height(ADV_INSVELX_dt_nums),width(ADV_INSVELX_dt_nums));

% Test 15b - Rate of Change YVEL
    INSVELY_dt_max = 1;
    INSVELY_dt_min = NaN;
    
    ADV_INSVELY_dt_nums = Vector_nums(:,2);

    Flagged_Struct.INSVELY_dt = ones(height(ADV_INSVELY_dt_nums),width(ADV_INSVELY_dt_nums));

 % repeat numbers
 %ADV_INSVELX_rpt_nums = Vector_nums(:,1);
% % Test 16a - X-Velocity Spike
    INSVELX_spk_max = 1;  % m/s
    INSVELX_spk_min = NaN;
    
    ADV_INSVELX_spk_nums = Vector_nums(:,1);

    Flagged_Struct.INSVELX_spk = ones(height(ADV_INSVELX_spk_nums),width(ADV_INSVELX_spk_nums));

% % Test 16b - Y-Velocity Spike
    INSVELY_spk_max = 1;
    INSVELY_spk_min = NaN;
    
    ADV_INSVELY_spk_nums = Vector_nums(:,2);

    Flagged_Struct.INSVELY_spk = ones(height(ADV_INSVELY_spk_nums),width(ADV_INSVELY_spk_nums));


%% Standard deviation filter for each velocity

Flagged_Struct.INSVEL_stdfltr = ones(height(ADV_INSVELY_spk_nums),width(ADV_INSVELY_spk_nums));

%% Calculate 2-D Current Speed and Direction


ADV_VEL =  Vector_nums(:,1:3);

[ADV_Current_Magnitude(:,1),ADV_Current_Direction(:,1),heading] = f_ADV_Speed_Direction_Pat(ADV_VEL);

% if PLOT == 1

%    figure
%        plot(Vector_datetime_UTC(:,1),ADV_Current_Direction(:,1),'.')

%    figure
%        plot(Vector_datetime_UTC(:,1),ADV_Current_Magnitude(:,1),'.')

 %        polarplot(deg2rad(ADV_Current_Direction(:,1)),ADV_Current_Magnitude(:,1),'.')
 % end
 
 %%   *******************     TESTS ******************************

 
%% Test 6
clear SIGSTR_Flagged

for i = 1:height(ADV_Amp1_nums)
%     if isnan(ADV_Amp1_nums(i,1)) || isnan(ADV_Amp2_nums(i,1)) || ...
%         isnan(ADV_Amp3_nums(i,1))  
%        Flagged_Struct.SIGSTR(i,1) = 0; % Fail
% 
%     else
        if ADV_Amp1_nums(i,1) >= SIGSTR_min && ADV_Amp1_nums(i,1) <= SIGSTR_max && ...
            ADV_Amp2_nums(i,1) >= SIGSTR_min && ADV_Amp2_nums(i,1) <= SIGSTR_max && ...
            ADV_Amp3_nums(i,1) >= SIGSTR_min && ADV_Amp3_nums(i,1) <= SIGSTR_max 

            Flagged_Struct.SIGSTR(i,1) = 1; % PASS
        else
            Flagged_Struct.SIGSTR(i,1) = 0; % FAIL
        end
%    end
end  

% Issue with the next line of code discrepancy 
%  Flagged_Struct.SIGSTR((ADV_Amp1_nums <= SIGSTR_min | ADV_Amp1_nums >= SIGSTR_max) & ...
%                        (ADV_Amp2_nums <= SIGSTR_min | ADV_Amp2_nums >= SIGSTR_max) & ...
%                        (ADV_Amp3_nums <= SIGSTR_min | ADV_Amp3_nums >= SIGSTR_max) ) = 0;
%             
            
 SIGSTR_Flagged = find(Flagged_Struct.SIGSTR(:,1)==0);    

%SIGSTR
if PLOT == 1

    str1 = strcat('SIGSTR all channels ',m);  % Method 1

    figure (5)
    plot(Vector_datetime_UTC(:,1),ADV_Amp1_nums(:,1),'b .')
    hold on
    plot(Vector_datetime_UTC(:,1),ADV_Amp2_nums(:,1),'r .')
    hold on
    plot(Vector_datetime_UTC(:,1),ADV_Amp3_nums(:,1),'g .')

        legend('Beam 1','Beam 2','Beam 3')
        title(str1)  
        ylabel(' ')
        xlabel('TIME UTC')
    hold on
    yline(SIGSTR_max)
    hold on
    yline(SIGSTR_min)
    
    figure (5)
    hold on
    plot(Vector_datetime_UTC(SIGSTR_Flagged,1),ADV_Amp1_nums(SIGSTR_Flagged,1),'r *')
    hold on
    plot(Vector_datetime_UTC(SIGSTR_Flagged,1),ADV_Amp2_nums(SIGSTR_Flagged,1),'r *')
    hold on
    plot(Vector_datetime_UTC(SIGSTR_Flagged,1),ADV_Amp3_nums(SIGSTR_Flagged,1),'r *')
         hold on
        legend('Beam 1','Beam 2','Beam 3','Flagged')    
end    
% Capture Points Removed Data
Test6_pts_flagged = length(SIGSTR_Flagged);

% Check if points have already been identified in previous tests 
[Total_pts_flagged_post_6,ADV_MasterFlag] = f_ADV_Master_Flag_Matrix(Flagged_Struct);


%% Test 7
clear Corr_Flagged

for i = 1:height(ADV_Corr1_nums)
%     if isnan(ADV_Corr1_nums(i,1)) || isnan(ADV_Corr2_nums(i,1)) || ...
%         isnan(ADV_Corr3_nums(i,1))  
% 
%        Flagged_StructCorr(i,1) = 0; % FAIL

%    else
        if ADV_Corr1_nums(i,1) >= Corr_min && ADV_Corr1_nums(i,1) <= Corr_max && ...
            ADV_Corr2_nums(i,1) >= Corr_min && ADV_Corr2_nums(i,1) <= Corr_max && ...
            ADV_Corr3_nums(i,1) >= Corr_min  && ADV_Corr3_nums(i,1) <= Corr_max 

            Flagged_Struct.Corr(i,1) = 1; % PASS
        else
            Flagged_Struct.Corr(i,1) = 0; % FAIL
        end
%    end

end    

% Issue with the next line of code discrepancy 
% Flagged_Struct.Corr( ADV_Corr1_nums <= Corr_min | ADV_Corr1_nums >= Corr_max & ...
%                      ADV_Corr2_nums <= Corr_min | ADV_Corr2_nums >= Corr_max & ...
%                      ADV_Corr3_nums <= Corr_min | ADV_Corr3_nums >= Corr_max) = 0; 

Corr_Flagged = find(Flagged_Struct.Corr(:,1)==0);    


str1 = strcat('Correlation all channels ', m);  % Method 1

% Correlation
if PLOT == 1

    figure (6)
    plot(Vector_datetime_UTC(:,1),ADV_Corr1_nums(:,1),'b .')
    hold on
    plot(Vector_datetime_UTC(:,1),ADV_Corr2_nums(:,1),'r .')
    hold on
    plot(Vector_datetime_UTC(:,1),ADV_Corr3_nums(:,1),'g .')

        legend('Beam 1','Beam 2','Beam 3')
        title(str1)  
        ylabel(' ')
        xlabel('TIME UTC')

    figure (6)
    hold on
    plot(Vector_datetime_UTC(Corr_Flagged,1),ADV_Corr1_nums(Corr_Flagged,1),'r .')
    %hold on
    %plot(Vector_datetime_UTC(Corr_Flagged,1),ADV_Corr2_nums(Corr_Flagged,1),'r .')
    %hold on
    %plot(Vector_datetime_UTC(Corr_Flagged,1),ADV_Corr3_nums(Corr_Flagged,1),'r .')
    %     hold on
        legend('Beam 1','Beam 2','Beam 3','Flagged')    
end

% Capture Points Removed Data
Test7_pts_flagged = length(Corr_Flagged);

% Check if points have already been identified in previous tests 
[Total_pts_flagged_post_7,ADV_MasterFlag] = f_ADV_Master_Flag_Matrix(Flagged_Struct);


 %% Test 10 
clear Current_Mag_Flagged

% for i = 1:height(ADV_Current_Magnitude)
%     for j = 1:width(ADV_Current_Magnitude)
%         if isnan(ADV_Current_Magnitude(i,j))
%             Flagged_Struct.Horz_Current_Mag(i,j) = 1; % PASS
%         else
%             if ADV_Current_Magnitude(i,j)  >= ADV_Current_Magnitude_min && ADV_Current_Magnitude(i,j) <= ADV_Current_Magnitude_max 
%                 Flagged_Struct.Horz_Current_Mag(i,j) = 1; % PASS
%             else
%                 Flagged_Struct.Horz_Current_Mag(i,j) = 0; % FAIL
%             end
%         end
%     end
% end    

Flagged_Struct.Horz_Current_Mag(ADV_Current_Magnitude <= ADV_Current_Magnitude_min | ADV_Current_Magnitude >= ADV_Current_Magnitude_max ) = 0;
 

Current_Mag_Flagged = find(Flagged_Struct.Horz_Current_Mag(:,1)==0);    


str1 = strcat('Horizontal Current Magnitude (Points removed from Magnitude criteria)', m);  % Method 1

if PLOT == 1

    figure (8)
        plot(Vector_datetime_UTC(:,1),ADV_Current_Magnitude(:,1),'b .')
        title(str1)  
        yline(0)
        ylabel('m/s ')
        xlabel('TIME UTC')
    figure(8)
        hold on
        plot(Vector_datetime_UTC(Current_Mag_Flagged,1),ADV_Current_Magnitude(Current_Mag_Flagged,1),'r .')
end
% Capture Points Removed Data
Test10_pts_flagged = length(Current_Mag_Flagged);

% Check if points have already been identified in previous tests 
[Total_pts_flagged_post_10,ADV_MasterFlag] = f_ADV_Master_Flag_Matrix(Flagged_Struct);
    
%% Test 11
clear Current_Dir_Flagged

% for i = 1:height(ADV_Current_Direction)
%     for j = 1:width(ADV_Current_Direction)
%         if isnan(ADV_Current_Direction(i,j))
%            Flagged_Struct.Horz_Current_Dir(i,j) = 1; % PASS
%         else
%             if ADV_Current_Direction(i,j)  >= ADV_Current_Direction_min && ADV_Current_Direction(i,j) <= ADV_Current_Direction_max 
%                 Flagged_Struct.Horz_Current_Dir(i,j) = 1; % PASS
%             else
%                Flagged_Struct.Horz_Current_Dir(i,j) = 0; % FAIL
%             end
%         end
%     end
% end    

% Flagged_Struct.Horz_Current_Dir(ADV_Current_Direction <= ADV_Current_Direction_min | ADV_Current_Direction >= ADV_Current_Direction_max) = 0; 
% 
% Current_Dir_Flagged = find(Flagged_Struct.Horz_Current_Dir(:,1)==0);    
% 
% 
% str1 = strcat('Horizontal Current Magnitude (Points removed from Directional criteria)', m);  % Method 1

if PLOT == 1

%     figure (9)
%         plot( Vector_datetime_UTC(:,1),ADV_Current_Direction(:,1),'b .')
%         title(str1)  
%         yline(0)
%         ylabel('m/s ')
%         xlabel('TIME UTC')
%     figure(9)
%         hold on
%         plot(Vector_datetime_UTC(Current_Dir_Flagged,1),ADV_Current_Direction(Current_Dir_Flagged,1),'r .')
%         
end
% Capture Points Removed Data
%Test11_pts_flagged = length(Current_Dir_Flagged);

% Check if points have already been identified in previous tests 
%[Total_pts_flagged_post_11,ADV_MasterFlag] = f_ADV_Master_Flag_Matrix(Flagged_Struct);

    
%% Test 12a Max XVEL
clear INSVELX_max_Flagged

% for i = 1:height(INSVELX_max_nums)
%     for j = 1:width(INSVELX_max_nums)
%         if isnan(INSVELX_max_nums(i,j))
%             Flagged_Struct.INSVELX_max(i,j) = 1; % PASS
%         else
%             if abs(INSVELX_max_nums(i,j)) <= INSVELX_max_max 
%                 Flagged_Struct.INSVELX_max(i,j) = 1; % PASS
%             else
%                 Flagged_Struct.INSVELX_max(i,j) = 0; % FAIL
%             end
%         end
%     end
% end    

Flagged_Struct.INSVELX_max(abs(INSVELX_max_nums) >= INSVELX_max_max) = 0;

INSVELX_max_Flagged = find(Flagged_Struct.INSVELX_max(:,1)==0);    

str1 = strcat('Test 12a INS VEL X Max', m);  % Method 1

if PLOT == 1

    figure(10)
        plot(Vector_datetime_UTC(:,1),Vector_nums(:,1), 'b .')
        title(str1)  
        yline(0)
        ylabel('m/s ')
        xlabel('TIME UTC')
    figure(10)
        hold on
        plot(Vector_datetime_UTC(INSVELX_max_Flagged,1),Vector_nums(INSVELX_max_Flagged,1),'r .')
end

% Capture Points Removed Data
Test12a_pts_flagged = length(INSVELX_max_Flagged);

% Check if points have already been identified in previous tests 
[Total_pts_flagged_post_12a,ADV_MasterFlag] = f_ADV_Master_Flag_Matrix(Flagged_Struct);


 %% Test 12b Max YVEL
 clear INSVELY_max_Flagged

% for i = 1:height(INSVELY_max_nums)
%     for j = 1:width(INSVELY_max_nums)
%         if isnan(INSVELY_max_nums(i,j))
%             Flagged_Struct.INSVELY_max(i,j) = 1; % PASS
%         else
%             if abs(INSVELY_max_nums(i,j)) <= INSVELY_max_max 
%                 Flagged_Struct.INSVELY_max(i,j) = 1; % PASS
%             else
%                 Flagged_Struct.INSVELY_max(i,j) = 0; % FAIL
%             end
%         end
%     end
% end   

Flagged_Struct.INSVELY_max(abs(INSVELY_max_nums) >= INSVELY_max_max) = 0;

INSVELY_max_Flagged = find(Flagged_Struct.INSVELY_max(:,1)==0);    

str1 = strcat('Test 12b INS VEL Y Max', m);  % Method 1

if PLOT == 1

    figure(11)
        plot(Vector_datetime_UTC(:,1),Vector_nums(:,2), 'b *')
        title(str1)  
        yline(0)
        ylabel('m/s ')
        xlabel('TIME UTC')
    figure(11)
        hold on
        plot(Vector_datetime_UTC(INSVELY_max_Flagged,1),Vector_nums(INSVELY_max_Flagged,2),'r .')
end

% Capture Points Removed Data
Test12b_pts_flagged = length(INSVELY_max_Flagged);

% Check if points have already been identified in previous tests 
[Total_pts_flagged_post_12b,ADV_MasterFlag] = f_ADV_Master_Flag_Matrix(Flagged_Struct);


 %% Test 13 Max ZVEL
% Vertical Velocity -  MIN=0, MAX-X=3000 [cm/s], 
clear INSVELZ_max_Flagged

% for i = 1:height(INSVELZ_max_nums)
%     for j = 1:width(INSVELZ_max_nums)
%         if isnan(INSVELZ_max_nums(i,j))
%             Flagged_Struct.INSVELZ_max(i,j) = 1; % PASS
%         else
    %             if abs(INSVELZ_max_nums(i,j)) <= INSVELZ_max_max 
    %                 Flagged_Struct.INSVELZ_max(i,j) = 1; % PASS
    %             else
    %                 Flagged_Struct.INSVELZ_max(i,j) = 0; % FAIL
    %             end
%         end
%     end
% end    

Flagged_Struct.INSVELZ_max(abs(INSVELZ_max_nums) >= INSVELZ_max_max) = 0;

INSVELZ_max_Flagged = find(Flagged_Struct.INSVELZ_max(:,1)==0);    

str1 = strcat('Test 13 INS VEL Z Max', m);  % Method 1

if PLOT == 1

    figure(12)
        plot(Vector_datetime_UTC(:,1),Vector_nums(:,3), 'b *')
        title(str1)  
        yline(0)
        ylabel('m/s ')
        xlabel('TIME UTC')
    figure(12)
        hold on
        plot(Vector_datetime_UTC(INSVELZ_max_Flagged,1),Vector_nums(INSVELZ_max_Flagged,3),'r +')
end

% Capture Points Removed Data
Test13_pts_flagged = length(INSVELZ_max_Flagged);
    
% Check if points have already been identified in previous tests 
[Total_pts_flagged_post_13,ADV_MasterFlag] = f_ADV_Master_Flag_Matrix(Flagged_Struct);


%%  Test 15a
% Rate of Change x velocity
% Test 15: u, v Rate of Change (for each bin)
% 	Instrument Coordinate Velocity (abs val of x-vel at t minus x-vel at t+1, same for y-vel), SUSPECT=800 [cm/s], FAIL=1200 [cm/s]
clear INSVELX_dt_Flagged 

for i = 1:height(ADV_INSVELX_dt_nums)-1
%     if isnan(ADV_INSVELX_dt_nums(i,1)) || isnan(ADV_INSVELX_dt_nums(i+1,1))  
%                 Flagged_Struct.INSVELX_dt(i,1) = 1; % PASS
%         else
                if abs(ADV_INSVELX_dt_nums(i,1) - ADV_INSVELX_dt_nums(i+1,1))  <= INSVELX_dt_max   
                    Flagged_Struct.INSVELX_dt(i,1) = 1; % PASS
                else
                    Flagged_Struct.INSVELX_dt(i,1) = 0; % FAIL
                end
%     end
 end    


INSVELX_dt_Flagged = find(Flagged_Struct.INSVELX_dt(:,1)==0);    

if PLOT == 1
    str1 = strcat('Test 15a INS VEL X Rate of Change ', m);  % Method 1

    figure (13)
        plot(Vector_datetime_UTC(:,1),Vector_nums(:,1),'b .')
        title(str1)  
        yline(0)
        ylabel('m/s ')
        xlabel('TIME UTC')
    figure(13)
        hold on
        plot(Vector_datetime_UTC(INSVELX_dt_Flagged,1),Vector_nums(INSVELX_dt_Flagged,1),'r .')
end
% Capture Points Removed Data
Test15a_pts_flagged = length(INSVELX_dt_Flagged);

% Check if points have already been identified in previous tests 
[Total_pts_flagged_post_15a,ADV_MasterFlag] = f_ADV_Master_Flag_Matrix(Flagged_Struct);

%% Test 15b
% Rate of Change y velocity
% Test 15: u, v Rate of Change (for each bin)
% 	Instrument Coordinate Velocity (abs val of x-vel at t minus x-vel at t+1, same for y-vel), SUSPECT=800 [cm/s], FAIL=1200 [cm/s]
clear  INSVELY_dt_Flagged 
 
for i = 1:height(ADV_INSVELY_dt_nums)-1
    
%    if isnan(ADV_INSVELY_dt_nums(i,1)) || isnan(ADV_INSVELY_dt_nums(i+1,1))  
% 
%         Flagged_Struct.INSVELY_dt(i,1) = 1; % PASS
%     else    
            if abs(ADV_INSVELY_dt_nums(i,1) - ADV_INSVELY_dt_nums(i+1,1))  <= INSVELY_dt_max   
                Flagged_Struct.INSVELY_dt(i,1) = 1; % PASS
            else
                Flagged_Struct.INSVELY_dt(i,1) = 0; % FAIL
            end
%    end
    
end    

INSVELY_dt_Flagged = find(Flagged_Struct.INSVELY_dt(:,1)==0);    

% yvel
if PLOT == 1
    str1 = strcat('Test 15b INS VEL Y Rate of Change', m );  % Method 1

    figure (14)
        plot(Vector_datetime_UTC(:,1),Vector_nums(:,2),'b .')
        title(str1)  
        yline(0)
        ylabel('m/s ')
        xlabel('TIME UTC')
    figure(14)
        hold on
        plot(Vector_datetime_UTC(INSVELY_dt_Flagged,1),Vector_nums(INSVELY_dt_Flagged,1),'r .')
end

% Capture Points Removed Data
Test15b_pts_flagged = length(INSVELY_dt_Flagged);

% Check if points have already been identified in previous tests 
[Total_pts_flagged_post_15b,ADV_MasterFlag] = f_ADV_Master_Flag_Matrix(Flagged_Struct);

%% Test 16a
% XVelocity Spike
% Test 16: u, v Spike (for each bin)
% 	Instrument Coordinate Velocity 
% (abs val of x-vel at t minus average of x-vel at t+1 and t-1, same for y-vel), SUSPECT=800 [cm/s], FAIL=1200 [cm/s]
clear INSVELX_spk_Flagged 

for i = 1:height(ADV_INSVELX_spk_nums)-1  % over time

    if i ==1 
        Flagged_Struct.INSVELX_spk(i,1) = 1; % skip 1st time point
    else
%         if isnan(ADV_INSVELX_spk_nums(i-1,1)) || isnan(ADV_INSVELX_spk_nums(i,1)) || isnan(ADV_INSVELX_spk_nums(i+1,1))
%             Flagged_Struct.INSVELX_spk(i,1) = 1; % skip points surronded by NaN's
%         else
                if (abs(ADV_INSVELX_spk_nums(i,1)) - abs(((ADV_INSVELX_spk_nums(i-1,1)+ ADV_INSVELX_spk_nums(i+1,1))/2))) <= INSVELX_spk_max

                    Flagged_Struct.INSVELX_spk(i,1) = 1; % PASS
                else
                    Flagged_Struct.INSVELX_spk(i,1) = 0; % FAIL
                end
%         end
     end

end    

INSVELX_spk_Flagged = find(Flagged_Struct.INSVELX_spk(:,1)==0);    

if PLOT == 1

    str1 = strcat('Test 16a INS VEL X Spike', m);  % Method 1

    figure (15)
        plot(Vector_datetime_UTC(:,1),Vector_nums(:,1),'b .')
        title(str1)  
        yline(0)
        ylabel('m/s ')
        xlabel('TIME UTC')    
    figure(15)
        hold on
        plot(Vector_datetime_UTC(INSVELX_spk_Flagged,1),Vector_nums(INSVELX_spk_Flagged,1),'r .')
end

% Capture Points Removed Data
Test16a_pts_flagged = length(INSVELX_spk_Flagged);

% Check if points have already been identified in previous tests 
[Total_pts_flagged_post_16a,ADV_MasterFlag] = f_ADV_Master_Flag_Matrix(Flagged_Struct);

%% Test 16b
% YVelocity Spike
% Test 16: u, v Spike (for each bin)
% 	Instrument Coordinate Velocity 
% (abs val of x-vel at t minus average of x-vel at t+1 and t-1, same for y-vel), SUSPECT=800 [cm/s], FAIL=1200 [cm/s]
clear INSVELY_spk_Flagged 

for i = 1:height(ADV_INSVELY_spk_nums)-1

    if i ==1 
        Flagged_Struct.INSVELY_spk(i,1) = 1; % skip 1st time point
    else
        if isnan(ADV_INSVELY_spk_nums(i-1,1)) || isnan(ADV_INSVELY_spk_nums(i,1)) || isnan(ADV_INSVELY_spk_nums(i+1,1))
            Flagged_Struct.INSVELY_spk(i,1) = 1; % skip points surronded by NaN's
        else
            if (abs(ADV_INSVELY_spk_nums(i,1)) - abs(((ADV_INSVELY_spk_nums(i-1,1)+ ADV_INSVELY_spk_nums(i+1,1))/2))) <= INSVELY_spk_max
                Flagged_Struct.INSVELY_spk(i,1) = 1; % PASS
            else
                Flagged_Struct.INSVELY_spk(i,1) = 0; % FAIL
            end
        end
    end

end    

INSVELY_spk_Flagged = find(Flagged_Struct.INSVELY_spk(:,1)==0);    

if PLOT == 1

    str1 = strcat('Test 16b INS VEL Y Spike', m);  % Method 1

    figure (16)
        plot(Vector_datetime_UTC(:,1),Vector_nums(:,2),'b .')
        title(str1)  
        yline(0)
        ylabel('m/s ')
        xlabel('TIME UTC')
    figure(16)
        hold on
        plot(Vector_datetime_UTC(INSVELY_spk_Flagged,1),Vector_nums(INSVELY_spk_Flagged,2),'r .')
end

% Capture Points Removed Data
Test16b_pts_flagged = length(INSVELY_spk_Flagged);

% Check if points have already been identified in previous tests 
[Total_pts_flagged_post_16b,ADV_MasterFlag] = f_ADV_Master_Flag_Matrix(Flagged_Struct);

%% Test 17 - Repeat Values
% 	Instrument Coordinate Velocity – X, SUSPECT=3, FAIL=5
 %computed on xvel only but will remove values from all velocity locaations when included in master flagged table                
% Test 17: Flat Line
% 	Instrument Coordinate Velocity – X, SUSPECT=3, FAIL=5
% clear INSVELX_rpt_Flagged
% 
% for i = 1:height(ADV_INSVELX_rpt_nums)-1
%     for j = 1:width(ADV_INSVELX_rpt_nums)
%         if isnan(ADV_INSVELX_rpt_nums(i,j)) || isnan(ADV_INSVELX_rpt_nums(i+1,j))
%             Flagged_Struct.INSVELX_rpt(i,j) = 1; % PASS
%         else    
%             if ADV_INSVELX_rpt_nums(i,j) ~= ADV_INSVELX_rpt_nums(i+1,j)   
%                 Flagged_Struct.INSVELX_rpt(i,j) = 1; % PASS
%             else
%                 Flagged_Struct.INSVELX_rpt(i,j) = 0; % FAIL
%             end
%         end
%     end
% end    
% 
% INSVELX_rpt_Flagged = find(Flagged_Struct.INSVELX_rpt(:,1)==0);    
% 
% str1 = strcat('Test 17 ADV Repeated Values', m);  % Method 1
% 
% figure (17)
%     plot(Vector_datetime_UTC(:,1),Vector_nums(:,1),'b .')
%     title(str1)  
%     yline(0)
%     ylabel('m/s ')
%     xlabel('TIME UTC')
% figure(17)
%     hold on
%     plot(Vector_datetime_UTC(INSVELX_rpt_Flagged,1),Vector_nums(INSVELX_rpt_Flagged,1),'r .')
% 
% % Capture Points Removed Data
% Test17_pts_flagged = length(INSVELX_rpt_Flagged);
% 
% % Check if points have already been identified in previous tests 
% [Total_pts_flagged_post_17,ADV_MasterFlag] = f_ADV_Master_Flag_Matrix(Flagged_Struct);
% 
% 
 
%% See how new velocities looks

%if PLOT == 1

    str1 = strcat('Prior to QC (XVEL shown)', m);  % Method 1
    figure

    hold on
        plot(Vector_datetime_UTC(:,1),Vector_nums(:,1), 'b .')

        title(str1)  
        yline(0)
        ylabel('m/s ')
        xlabel('TIME UTC')
        ylim([-3 3])

%end


    Vector_nums_new = ones(length(Vector_nums),width(Vector_nums));

for i = 1:length(Vector_nums)

    if ADV_MasterFlag(i) == 1

        Vector_nums_new(i,:) = Vector_nums(i,:);

    else

         Vector_nums_new(i,:) = NaN;
    end

end



if PLOT == 1

    figure

    str1 = strcat('Results After QC (XVEL shown)', m);  % Method 1

    hold on
        plot(Vector_datetime_UTC(:,1),Vector_nums_new(:,1), 'b .')
        title(str1)  
        yline(0)
        ylabel('m/s ')
        xlabel('TIME UTC')
        ylim([-3 3])
end

%% Perfom Basic standard deviation filter 3 times for each velocity 

% ntimes     = 3;
% 
% %X Velocity
% ts  = Vector_nums_new(:,1);
% 
% %Remove the mean from the veocity ts
% 
% N = length(ts);
% [ts_mean,ts_var,ts_stdev] = f_stats_m_var_stdev(ts,N);
% 
% if ts_mean<0
% 
%     ts = ts+abs(ts_mean);   
% elseif ts_mean>0
%     ts = ts-ts_mean;   
% end
% 
% %Apply number of standard deviations for cutoff based on which ADV was used
% 
%       if k == k1 
%             no = 4.5;
%       else 
%             no = 3;
%       end
% 
% %Compute filtered ts
% window = sf*60*4;
% [xvel_ts_filter,xvel_mean,xvel_var,xvel_stdev,xvel_Ngood] = f_stdev_filter(ts,no,ntimes);
% [xvel_ts_filter] = f_moving_std_filter(Vector_nums_new(:,1),sf,window,3);
% xvel_ts_filter = xvel_ts_filter';
% n = num2str(no);
% 
% if PLOT == 1
% 
%     str1 = strcat(m,' Results After QC and std-dev filter ran 3 times ');  % Method 1
%     str2 = strcat(' (XVEL shown)',' number of std devs: ',n);
%     str3 = {str1,str2};
% 
%     figure
% 
%         plot(Vector_datetime_UTC(:,1),ts(:,1), 'b .')
%         hold on
%         plot(Vector_datetime_UTC(:,1),xvel_ts_filter, 'r .')
%         title(str3)  
%         yline(0)
%         ylabel('m/s ')
%         xlabel('TIME UTC')
% end
% 
% %Y Velocity    
% ts         = Vector_nums_new(:,2);
% 
% %Remove the mean from the veocity ts
% 
% N = length(ts);
% [ts_mean,ts_var,ts_stdev] = f_stats_m_var_stdev(ts,N);
% 
% if ts_mean<0
% 
%     ts = ts+abs(ts_mean);   
% elseif ts_mean>0
%     ts = ts-ts_mean;   
% end
% 
% %Apply number of standard deviations for cutoff based on which ADV was used
% 
%       if k == k1 
%             no = 11;
%       else 
%             no = 4.5;
%       end
% 
% %Compute filtered ts
% 
% [yvel_ts_filter,yvel_mean,yvel_var,yvel_stdev,yvel_Ngood] = f_stdev_filter(ts,3,ntimes);
% [yvel_ts_filter] = f_moving_std_filter(Vector_nums_new(:,2),sf,window,3);
% yvel_ts_filter = yvel_ts_filter';
% if PLOT == 1
%     n = num2str(no);
% 
%     str1 = strcat(m,' Results After QC and std-dev filter ran 3 times');  % Method 1
%     str2 = strcat(' (YVEL shown)',' std deviations: ',n);
%     str3 = {str1,str2};
% 
%     figure
%         hold on
%         plot(Vector_datetime_UTC(:,1),ts(:,1), 'b .')
%         plot(Vector_datetime_UTC(:,1),yvel_ts_filter, 'r .')
%         title(str3)  
%         yline(0)
%         ylabel('m/s ')
%         xlabel('TIME UTC')
% 
% end
% %Z Velocity    
% ts         = Vector_nums_new(:,3);
% 
% %Remove the mean from the veocity ts
% 
% N = length(ts);
% [ts_mean,ts_var,ts_stdev] = f_stats_m_var_stdev(ts,N);
% 
% if ts_mean<0
% 
%     ts = ts+abs(ts_mean);   
% elseif ts_mean>0
%     ts = ts-ts_mean;   
% end
% 
% %Apply number of standard deviations for cutoff
% 
%       if k == k1 
%             no = 12;
%       else 
%             no = 9;
%       end
% 
% %Compute filtered ts
% 
% [zvel_ts_filter,zvel_mean,zvel_var,zvel_stdev,zvel_Ngood] = f_stdev_filter(ts,3,ntimes);
% [zvel_ts_filter] = f_moving_std_filter(Vector_nums_new(:,3),sf,window,3);
% zvel_ts_filter = zvel_ts_filter';
% 
% if PLOT == 1
% 
%     n = num2str(no);
% 
%     str1 = strcat( m,' Results After QC and standard deviation filter ran 3 times,');  % Method 1
%     str2 = strcat(' (ZVEL shown)',' std deviations: ',n);
%     str3 = {str1,str2};
%     figure
%         hold on
%         plot(Vector_datetime_UTC(:,1),ts(:,1), 'b .')
%         plot(Vector_datetime_UTC(:,1),zvel_ts_filter, 'r .')
%         title(str3)  
%         yline(0)
%         ylabel('m/s ')
%         xlabel('TIME UTC')
% end
% 
% %% Include in Master Flagged Vector
% 
% for i = 1:length(xvel_ts_filter)
% 
%     if isnan(xvel_ts_filter(i)) ||  isnan(yvel_ts_filter(i)) ||  isnan(zvel_ts_filter(i))
%         Flagged_Struct.INSVEL_stdfltr(i) = 0; % Pass
%     else
%         Flagged_Struct.INSVEL_stdfltr(i) = 1; % FAIL
%     end
% 
% end    
% 
% INSVEL_stdfltr_Flagged = find(Flagged_Struct.INSVEL_stdfltr(:,1)==0);  
% 
% if PLOT == 1
% 
%     str1 = strcat(m,' Test 21 ADV Std filter to each vel (XVEL shown)');  % Method 1
% 
%     figure 
%     hold on
%         plot(Vector_datetime_UTC(:,1),Vector_nums_new(:,1),'b .')
%         title(str1)  
%         yline(0)
%         ylabel('m/s ')
%         xlabel('TIME UTC')
%     hold on
%     plot(Vector_datetime_UTC(INSVEL_stdfltr_Flagged,1),Vector_nums_new(INSVEL_stdfltr_Flagged,1),'r .')
% 
% end

% Capture Points Removed Data
Test21_pts_flagged = 0; %length(INSVEL_stdfltr_Flagged);

% Check if points have already been identified in previous tests 
[Total_pts_flagged_post_21,ADV_MasterFlag] = f_ADV_Master_Flag_Matrix(Flagged_Struct);
Total_pts_flagged_post_21=0;
Test21_pts_flagged=0;

%% Build a table to show results

% delta5 =  Total_pts_flagged_post_7 - Total_pts_flagged_post_6;
% delta6 =  Total_pts_flagged_post_10 - Total_pts_flagged_post_7;
% delta7 =  Total_pts_flagged_post_11 - Total_pts_flagged_post_10;
% delta8 = Total_pts_flagged_post_12a - Total_pts_flagged_post_11;
% delta9 = Total_pts_flagged_post_12b - Total_pts_flagged_post_12a;
% delta10 = Total_pts_flagged_post_13 - Total_pts_flagged_post_12b;
% delta11 = Total_pts_flagged_post_15a - Total_pts_flagged_post_13;
% delta12 =  Total_pts_flagged_post_15b - Total_pts_flagged_post_15a;
% delta13 =  Total_pts_flagged_post_16a - Total_pts_flagged_post_15b;
% delta14 =  Total_pts_flagged_post_16b - Total_pts_flagged_post_16a;
% delta15 =  Total_pts_flagged_post_21 - Total_pts_flagged_post_16b;

delta5 =  Total_pts_flagged_post_7 - Total_pts_flagged_post_6;
delta6 =  Total_pts_flagged_post_10 - Total_pts_flagged_post_7;
delta7 =  Total_pts_flagged_post_12a - Total_pts_flagged_post_10;
%delta8 = Total_pts_flagged_post_12a - Total_pts_flagged_post_11;
delta9 = Total_pts_flagged_post_12b - Total_pts_flagged_post_12a;
delta10 = Total_pts_flagged_post_13 - Total_pts_flagged_post_12b;
delta11 = Total_pts_flagged_post_15a - Total_pts_flagged_post_13;
delta12 =  Total_pts_flagged_post_15b - Total_pts_flagged_post_15a;
delta13 =  Total_pts_flagged_post_16a - Total_pts_flagged_post_15b;
delta14 =  Total_pts_flagged_post_16b - Total_pts_flagged_post_16a;
delta15 =  Total_pts_flagged_post_21 - Total_pts_flagged_post_16b;

Test6{1} = 'Signal Strength (db)';
Test6(1,(2:6)) = num2cell([ SIGSTR_min,  SIGSTR_max,  Test6_pts_flagged , Total_pts_flagged_post_6, (Total_pts_flagged_post_6/height(Vector_nums))*100 ]);

Test7{1} = 'Correlation (%)';
Test7(1,(2:6)) = num2cell([ Corr_min,  Corr_max,  Test7_pts_flagged , delta5,(Total_pts_flagged_post_7/height(Vector_nums))*100 ]);

Test10{1} = 'Current Magnitude (m/s)';
Test10(1,(2:6)) = num2cell([ ADV_Current_Magnitude_min,  ADV_Current_Magnitude_max,  Test10_pts_flagged , delta6, (Total_pts_flagged_post_10/height(Vector_nums))*100  ]);
%Test11{1} = 'Current Direction (deg)';
%Test11(1,(2:6)) = num2cell([ ADV_Current_Direction_min,  ADV_Current_Direction_max,  Test11_pts_flagged , delta7, (Total_pts_flagged_post_11/height(Vector_nums))*100 ]);

Test12a{1} = 'INS XVEL MAX (m/s)';
Test12a(1,(2:6)) = num2cell([INSVELX_max_min, INSVELX_max_max,  Test12a_pts_flagged , delta7, (Total_pts_flagged_post_12a/height(Vector_nums))*100 ]);

Test12b{1} = 'INS YVEL MAX (m/s)';
Test12b(1,(2:6)) = num2cell([INSVELY_max_min, INSVELY_max_max,  Test12b_pts_flagged , delta9, (Total_pts_flagged_post_12b/height(Vector_nums))*100]);

Test13{1} = 'INS ZVEL MAX (m/s)';
Test13(1,(2:6)) = num2cell([INSVELZ_max_min, INSVELZ_max_max,  Test13_pts_flagged , delta10, (Total_pts_flagged_post_13/height(Vector_nums))*100 ]);

Test15a{1} = 'Ins X Vel Rate of Change (m/s)';
Test15a(1,(2:6)) = num2cell([ INSVELX_dt_min,  INSVELX_dt_max,  Test15a_pts_flagged , delta11, (Total_pts_flagged_post_15a/height(Vector_nums))*100 ]);

Test15b{1} = 'Ins Y Vel Rate of Change (m/s)';
Test15b(1,(2:6)) = num2cell([ INSVELY_dt_min,  INSVELY_dt_max,  Test15b_pts_flagged , delta12, (Total_pts_flagged_post_15b/height(Vector_nums))*100]);

Test16a{1} = 'Ins X Vel Spike (m/s)';
Test16a(1,(2:6)) = num2cell([ INSVELX_spk_min,  INSVELX_spk_max,  Test16a_pts_flagged , delta13, (Total_pts_flagged_post_16a/height(Vector_nums))*100 ]);

Test16b{1} = 'Ins Y Vel Spike (m/s)';
Test16b(1,(2:6)) = num2cell([ INSVELY_spk_min,  INSVELY_spk_max,  Test16b_pts_flagged , delta14,(Total_pts_flagged_post_16b/height(Vector_nums))*100  ]);

if k == k2 
    Test21{1} = 'Standard Deviation Filter x=3std | y=4.5std | z=9std';
elseif k == k1 
    Test21{1} = 'Standard Deviation Filter x=4.5std | y=11std | z=12std';
end


Test21(1,(2:6)) = num2cell([NaN ,NaN,  Test21_pts_flagged , delta15,(Total_pts_flagged_post_21/height(Vector_nums))*100 ]);


T = [Test6; Test7; Test10; Test12a; Test12b; Test13; Test15a; Test15b; Test16a; Test16b; Test21]; %; Test11; Test15a; Test15b; Test16a; Test16b; Test17; Test20; Test21]; 
 

str1 = strcat('# of points flagged   ');  % Method 1
str2 = strcat('Increase in Total points Flagged ');  % Method 1

f = uifigure;
colnames = {  'Test Description', 'Lower Limit', 'Upper Limit',str1 ,str2,'% of Total Points Flmpstagged After each test' };       
H = uitable(f, 'Data', T, 'ColumnName', colnames,'ColumnWidth',{110},'Position',[20 75 500 225]);
s = uistyle('HorizontalAlignment','right');
addStyle(H,s,'column',1);
H.RowName = {'Test6','Test7','Test10', 'Test12a', 'Test12b', 'Test13','Test15a','Test15b','Test16a','Test16b','Test21'};%,,'Test11','Test17','Test20','Test21'};

%%

    Vector_nums_new = ones(length(Vector_nums),width(Vector_nums));
    
for i = 1:length(Vector_nums)
    
    if ADV_MasterFlag(i) == 1
        
        Vector_nums_new(i,:) = Vector_nums(i,:);
        
    else
        
         Vector_nums_new(i,:) = NaN;
    end
    
end

%if PLOT == 1

    figure
    % percent_removed = num2str(100*Total_pts_flagged_post_21/length(Vector{:,1}));
    percent_removed = num2str(100*Total_pts_flagged_post_21/length(Vector{:,2}));
    str1 = strcat( 'ADV ',str,' Results Before and After QC (XVEL shown)')%,percent_removed,'% of the data removed)');  % Method 1
    %ax1 = subplot(211);
    hold on
        %plot(Vector_datetime_UTC(:,1),Vector{:,2},'b.')
        plot(Vector_datetime_UTC(:,1),Vector_nums_new(:,1), 'b .');
        title(str1)  
        yline(0)
        ylabel('m/s ')
        xlabel('TIME UTC')
        
     %ax2 = subplot(212);
     figure
     plot(Vector_nums_new(:,1), 'b .');
     %linkaxes([ax1 ax2],'x')
    
%end

toc
%% Save the Results
tic
Vector_new(:,1) = array2table(Vector_datetime_UTC);
Vector_new(:,2:13) = array2table(Vector_nums_new(:,1:end));

Vector_new.Properties.VariableNames = Vector.Properties.VariableNames;
Vector_new.Properties.VariableUnits = Vector.Properties.VariableUnits;

% Vector_System_Data_new(:,1) = array2table(Vector_System_Data_datetime_UTC);
% Vector_System_Data_new(:,2:13) = Vector_System_Data(:,2:end);
% Vector_System_Data_new.Properties.VariableNames = Vector_System_Data.Properties.VariableNames;
% Vector_System_Data_new.Properties.VariableUnits = Vector_System_Data.Properties.VariableUnits;
% 
ADV_QC_Data.Velocity = Vector_new;
% ADV_QC_Data.SystemData = Vector_System_Data_new;
% 
%mkdir("../../ADV-Nortek Vector/8.17.23 Deployment/ADV"+str+"/QC Data")
%save("../../ADV-Nortek Vector/8.17.23 Deployment/ADV"+str+"/QC Data/Vector_"+str+"_UTC_PostQC_mod.mat",'-struct','ADV_QC_Data','-v7.3')
% "../../ADV-Nortek Vector/8.17.23 Deployment/ADV"+str+"/Converted Data/Vector_"+str+"_UTC.mat"
%save('../../Data_Files/MODAQ/Vector_Velocity_Data_1/Processed_Data/2021_11_04_to_12_15_Vector_1_UTC_PostQC.mat','-struct','ADV_QC_Data','-v7.3')

save('Vector_UTC_PostQC.mat','-struct','ADV_QC_Data','-v7.3')
toc


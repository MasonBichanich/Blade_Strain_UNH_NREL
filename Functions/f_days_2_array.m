function [Output_Data_table] = f_days_2_array(monthrange,dayrange,file_path_dailyfiles,file_name_format,current_filepath)
%[Output_Data_table] = untitled(monthrange,dayrange,file_path_dailyfiles,file_name_format)
%The purpose of this function is to load in DailyQC'd matlab plots over
%a select date range
% Input
%
% monthrange               =        ex. [11 11]
% dayrange                 =        ex. [11 11]
% file_path_dailyfiles     =        ex. '../../Data_Files/MODAQ/Vector_Velocity_Data_2/Processed_Data/Daily_Files_Post_QC'
% file_name_format         =        ex. '*_ADV2_Bow_Velocity_Data_Post_QC_UTC_NREL_Final.mat'
% current_filepath         =        ex. 'C:\Users\pwo1002\OneDrive - USNH\BOX\Patrick OByrne\Data_Analysis\Scripts'
%
% Outputs
% 
% Output_Data_table        = combined table array over input range

%% load in QC'd Data over specified date range

% select the range of months and days you want
%monthrange = [11 11];
%dayrange = [09 24];

% change directory to Processed Bow ADV Data
%file_path_dailyfiles = '../../Data_Files/MODAQ/Vector_Velocity_Data_2/Processed_Data/Daily_Files_Post_QC';
cd(file_path_dailyfiles)

%file_name_format = '*_ADV2_Bow_Velocity_Data_Post_QC_UTC_NREL_Final.mat';
%%
unsorted_data = dir(file_name_format);
table = struct2table(unsorted_data);
sortedtable = sortrows(table, 'datenum');
sortedtable = sortedtable(1:end-1,:); % removes last day of data because the filename is not the same length as all the rest
file_dirmaster = table2struct(sortedtable);

filenames = vertcat(file_dirmaster.name); %lists file names in vertical vector

year_ = string(filenames(:,1:4)); %indexes positions of hour labels in dir list
month_ = string(filenames(:,6:7));
day_ = string(filenames(:,9:10));

% convert from string to number
day_num = str2double(day_);
month_num = str2double(month_);

month_name = datestr(datetime(1,month_num,1),'mmm');

% Determine what indices in the month_ day_ year_ vectors fall within this
% range

j=1;
for i = 1:length(month_)
    % If days span only one month
    if ( month_num(i)== monthrange(1) &&  day_num(i)>=dayrange(1) && day_num(i)<=dayrange(end))
        good_indicies(j) = i;
        j = j+1;
    else
        i=i+1;
    end 
end
        
% Now need to load files that match these indices
j =1;
for i = 1:length(file_dirmaster)
    
    if (i>=good_indicies(1) && i<= good_indicies(end))
        labels(j,1) = strcat(month_name(i,:),'_',day_(i),'_',year_(i));
        Data(1).(labels(j,1)) = load(filenames(i,:));
        j=j+1;
    else
        i=i+1;
    end
end

% Combine all days into one data array
Output_Data_table = [];
for i =1:j-1
    Output_Data_table = vertcat(Output_Data_table,Data(1).(labels(i)).data(:,:));
end


% change directory back to script directory
cd(current_filepath)
end


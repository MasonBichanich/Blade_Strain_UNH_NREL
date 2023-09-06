function [phase,ave_data] = f_phase_averaging(data,time,dg,peak_ind)
% This function takes wave data (tested with sinusoidal wave) and creates a
% single, phase-averaged waveform. The function uses the peak as an
% arbitrary phase location marker, rather than the zero crosing (no need to
% remove low frequency trends. The parameters of the peakfinder function
% may need to be adjusted on a case-by-case basis

%  OUTPUTS
%     phase = vector of 0-360, length depends on averaging window chosen
%     ave_data = single waveform of averaged data
%
%  INPUTS
%     data = wave data that you'd like to phase average. This fuction will
%     tolerate NaN
%     time = time vector for the data. It needs to be numbers.
%     dg = number of degrees you like to average over. 1 degree is the
%     minimum
%     peak_ind = indicies of the data vector where you'd like to break up
%     the data
%break up the data at the peaks
broken = cell(length(peak_ind),1);
broken_time = cell(length(peak_ind),1);
broken{1} = data(1:peak_ind(1));
broken_time{1} = time(1:peak_ind(1));
for i = 2:length(peak_ind)
       broken{i} = data(peak_ind(i-1):peak_ind(i));
       broken_time{i} = time(peak_ind(i-1):peak_ind(i));
end

% throw away first and last bc they may be incomplete phases
broken(1) = [];
broken(end) = [];
broken_time(1) = [];
broken_time(end) = [];
% % check if it worked
figure; hold on;
for i = 1:length(broken)
    plot(broken_time{i},broken{i})
end
title('Data broken into phases based on peak location')
%% now lets average these mf phases. They are all different lengths. First try will be to make them all the same length 
% and Nan the missing data
% create a new compass, 0-360 and the length of each cell
% and use outerjoin

% making the compass, adding to the broken up data in the correct cell
for i = 1:length(broken)
    comp{i,:} = round(linspace(0,359,length(broken{i})));
    broken_comp{i,:} = [comp{i,:}',broken{i,:}];

end
%combining them all into a table (each degree has a bunch of nans and data
%points to average over
joined = array2table(broken_comp{1});
for i = 2:length(broken_comp)
    broken_comp{i} = array2table(broken_comp{i});
    joined = outerjoin(broken_comp{i},joined,Keys=1,MergeKeys=true);
end
joined = table2array(joined);
%%
% the data should be ensemble averaged BEFORE being phase averaged. Compass
% resolution is NOT 1 degree
ens_ave=[];
for i = 1:length(joined)
    for j = 1:length(joined(:,1))/dg
        ens_ave{i}(j,:) = mean(joined(dg*(j-1)+1:dg*j,i),'omitnan');
    end
end
ens_ave = cell2mat(ens_ave);
phase = ens_ave(:,1);
phase_averaged = mean(ens_ave(:,2:end),2,'omitnan');
pa_end = phase_averaged(round(3/4*length(phase_averaged)):end);
pa_start = phase_averaged(1:round(3/4*length(phase_averaged)-1));
ave_data = [pa_end;pa_start];
end
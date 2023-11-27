function [ts_filter,mean,var,stdev,Ngood] = f_stdev_filter(ts,no,ntimes,window)
%function [ts_filter,mean,var,stdev,Ngood] = f_stdev_filter(ts,no)
%This function deterimines if as value in a time series fall within within "no" # of
%standard deviations from the mean. This filter is repeated "ntimes" with
%the new mean and standard deviation of the iteratively filtered results
% Inputs
%      ts         = time series
%      no         = number of standard deviations
%      ntimes     = how many times you want to run the filter on your ts
% Outputs
%      ts_filter  = filtered time series
%      mean       = mean value of the filtered time series
%      Ngood      = Number of good values left in time series

N = length(ts);

[mean,var,stdev] = f_stats_m_var_stdev(ts);

for j =1:ntimes
    Ngood = 0;
    for i=1:N
    
    if ts(i) < mean + no*stdev  && ts(i) > mean - no*stdev  
        
        ts_filter(i,1) = ts(i);
        Ngood = Ngood +1;
    else 
        
        ts_filter(i,1) = NaN;
        
    end
    
end

[mean,var,stdev] = f_stats_m_var_stdev(ts_filter);

end

end


function [mean_val,var,stdev,good_points] = f_stats_m_var_stdev(ts,N)
%function [mean,var,stdev] = f_stats_m_var_stdev(ts,Ngood)
% This function calculates the mean, variance and standard deviation of the
% time series (ts) with Ngood number of data points ~= NaN
%Inputs
%   ts    = time series values
 %Outputs
%   mean = returns the average value of the time series
%   var  = variance of the time series
%   stdev = standard deviation of time series

N  = length(ts);

meansum=0;
good_points =1;

for i = 1:N 
    if ~isnan(ts(i))
    meansum = meansum + ts(i);
    good_points = good_points +1;
    else
        i = i+1;
    end
end

mean_val = meansum/good_points;

varsum = 0;
for i = 1:N 
    if ~isnan(ts(i))
    varsum = varsum + (ts(i)^2 - mean_val^2);
    else
        i = i+1;
    end
end

var = varsum/good_points;
stdev = sqrt(var);

 
end


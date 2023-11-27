function [filtered_ts] = f_moving_std_filter(ts,sf,new_T,no_std)
%This acts as a standard deviation filter. It evaluates the statistics in
%each window around each point, and uses those as filter parameters.
%
%ts = full time series (your data)
%t = your time data
%sf = sampling frequency
%new_T = the length of window you want to get the stats from (s)
%how many stdev away from mean you want to filter
%filtered_ts & filtered_t will be slightly shorter than original
%
%This function does NOT save the statistics calculated at each point

N=length(ts);
T=N/sf;
new_N=new_T*sf;
j=1;
% for i=new_N/2:1:N-new_N/2
%     new_ts=ts(i+1-new_N/2:i+new_N/2);
%     [mu(j,1),std(j,1)]=f_stats_m_var_stdev(new_ts,new_N);
%     j=j+1;
% end
mu=movmean(ts,new_N,'omitnan');
std=movstd(ts,new_N,'omitnan');
%% the filter
%moving standard deviation filter

for i=1:length(mu)
    if mu(i)-std(i)*no_std < ts(i) && ts(i) < mu(i)+std(i)*no_std
        filtered_ts(i)=ts(i);
       
    else
        filtered_ts(i)=NaN;

    end
end



end
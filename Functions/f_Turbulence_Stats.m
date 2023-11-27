function [mean_vel,anomoly,stress,TI,TKE] = Turbulence_Stats(uvw,sf,wind)
% This fuction pulls in your velocity data from an ADV and calculates
% various turbulence characteristics
% input variables:
% uvw   == matrix of velocities in the x,y,z directions (Nx3)  [m/s]
% t     == time   [s]
% sf    == sampling frequency   [Hz]
% wind  == window to average over (mean-removed stats)   [s]
% 

%% this is a sliding average. Should it be ensemble averages?

%% means and anomolies
for i = 1:3
    mean_vel(:,i) = movmean(uvw(:,i),wind*sf,'omitnan');
    anomoly(:,i) = uvw(:,i) - mean_vel(:,1);
end

%% Reynolds Stress

for i=1:6
    if i < 4
        stress{i} = movmean(anomoly(:,i).^2,wind*sf,'omitnan');
    elseif i == 4 || i == 5
        stress{i} = movmean(anomoly(:,1).*anomoly(:,i-2),sf*wind,'omitnan');
    else
        stress{i} = movmean(anomoly(:,2).*anomoly(:,3),sf*wind,'omitnan');
    end
end

%% Turbulence intensity

TI = 100*abs(sqrt(stress{1}(abs(mean_vel(:,1))>0.1))./mean_vel(abs(mean_vel(:,1))>0.1,1));

%% Turbulence Kinetic Energy

TKE = 0.5*(stress{1} + stress{2} + stress{3});

end
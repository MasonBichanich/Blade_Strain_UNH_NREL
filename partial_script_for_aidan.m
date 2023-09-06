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
%% making the figure with all the sensors
figure
for i = 1:8

    phase_averaged{i} = [ens_ave{i}(:,1),mean(ens_ave{i}(:,2:end),2,'omitnan')];

    %subtract the mean
    means{i} = mean(phase_averaged{i}(:,2));
    pa_a{i} = phase_averaged{i}(:,2)-means{i};

    % shifting the zero crossing since its arbitrary
    pa_end = pa_a{i}(round(1/4.5*length(phase_averaged{i})):end); 
    pa_start = pa_a{i}(1:round(1/4.5*length(phase_averaged{i})-1));
    pa_shifted{i} = [pa_end;pa_start];
    
    hold on
    plot(phase_averaged{i}(:,1),pa_shifted{i})
end
legend('S1','S2','S3','S4','S5','S6','S7','S8')
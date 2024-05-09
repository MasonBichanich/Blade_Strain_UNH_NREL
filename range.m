load('broken.mat');
load('broken_time');

%%
amp = {};
for i = 1:length(broken)
    for j = 1:length(broken{i})
         amp{i}(j) = table2array(abs(max(broken{i}{j})-min(broken{i}{j})))/2;
    end
end
amp_time={};
for i=1:length(broken_time)
    amp_time{i}=mean(broken_time{i});
end
%%

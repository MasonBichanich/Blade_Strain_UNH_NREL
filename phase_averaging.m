% load in strain data
A = exist('vlink_s1');

if A == 0
    vlink_8145 = readtable('../../Data_Files/MODAQ/Vlink/Processed_Data/8145_vLink_strain-2022-11-22.csv');
    vlink_28175 = readtable('../../Data_Files/MODAQ/Vlink/Processed_Data/28175_vLink_strain-2022-11-22.csv');

    upper_bound = 335729;  % this is where the turbine stops rotating
    % Extract vlink 8145 data
    vlink_8145 = vlink_8145(2:upper_bound,:);
    vlink_s1_s4_time = table2array(vlink_8145(:,1));
    vlink_s1_datenum = datenum(vlink_s1_s4_time);
    vlink_s1 = table2array(vlink_8145(:,2));
    vlink_s2 = table2array(vlink_8145(:,3));
    vlink_s3 = table2array(vlink_8145(:,4));
    vlink_s4 = table2array(vlink_8145(:,5));

    % Extract vlink 28175 data
    vlink_28175(1:2,:) = [];
    vlink_28175 = vlink_28175(1:upper_bound-1,:);
    vlink_s5_s8_time = table2array(vlink_28175(:,1));
    vlink_s5_datenum = datenum(vlink_s5_s8_time);
    vlink_s5 = table2array(vlink_28175(:,2));
    vlink_s6 = table2array(vlink_28175(:,3));
    vlink_s7 = table2array(vlink_28175(:,4));
    vlink_s8 = table2array(vlink_28175(:,5));
end

%% examine the data
plot(vlink_s1_datenum,vlink_s3)
title('Raw Data at S3')
zero = mean(vlink_s3);
hold on
yline(zero,'k')

%% chop up phases at peaks instead of at zero. The zero appears to have fluctuations over time, so I would need to do moving mean or something
% I've decided that the troughs of S3 are the most defined so I'll be using
% those
[dummy, peak_ind] = findpeaks(-vlink_s3,MinPeakDistance=100,MinPeakProminence=10);         %find the troughs as they are best defined
plot(vlink_s1_datenum(peak_ind),vlink_s3(peak_ind),'*r')                                     %confirm the peaks
title('Raw data truncated to turbine operation and highlighted peaks')

%% I'm going to see what happens if i remove the means from the data first
for i = 1:8
    if i < 5

        vlink_8145{:,i+1} = vlink_8145{:,i+1}-mean(vlink_8145{:,i+1},'omitnan');
    else
        vlink_28175{:,i-3} = vlink_28175{:,i-3}-mean(vlink_28175{:,i-3},'omitnan');
    end
end
%% break up the data at the troughs

%preallocation
broken = {cell(length(peak_ind),1), cell(length(peak_ind),1), cell(length(peak_ind),1), cell(length(peak_ind),1),...
    cell(length(peak_ind),1), cell(length(peak_ind),1), cell(length(peak_ind),1), cell(length(peak_ind),1),};
broken_time = cell(length(peak_ind),1);

%initializing
for i = 1:8
    if i < 5
        broken{i}{1} = vlink_28175(1:peak_ind(1),i+1);
    else
        broken{i}{1} = vlink_8145(1:peak_ind(1),i-3);
    end
end
broken_time{1} = vlink_s1_datenum(1:peak_ind(1));

%breaking the time
for i = 2:length(peak_ind)
       broken_time{i} = vlink_s1_datenum(peak_ind(i-1):peak_ind(i));
end

%breaking the data
for j = 1:8
    for i = 2:length(peak_ind)
        if j < 5
            broken{j}{i} = vlink_28175(peak_ind(i-1):peak_ind(i),j+1);
        else
            broken{j}{i} = vlink_8145(peak_ind(i-1):peak_ind(i),j-3);
        end
    end
end


% check if it worked
% for j = 1:8
%     figure
%     for i = 1:length(broken{1})
%         plot(broken_time{i},table2array(broken{j}{i}))
%         hold on
%         title(num2str(j))
%     end
% end

%% now lets average these mf phases. They are all different lengths. First try will be to make them all the same length 
% and Nan the missing data
% create a new compass, 0-360 and the length of each cell
% and use outerjoin

%create imaginary compass
for i = 1:height(broken{1})
    comp{i,:} = array2table(round(linspace(0,359,length(table2array(broken{1}{i}))))');
end

%combine imaginary compass with broken data
broken_comp = broken;
for j = 1:8
    for i = 1:height(broken{1})
        broken_comp{j}{i} = [comp{i},broken{j}{i}];
    end
end


%combining them all into a table (each degree has a bunch of nans and data
%points to average over) (this takes a long time)
tic
for j = 1:8
    joined{j} = (broken_comp{1}{1});
    for i = 2:height(broken_comp{1})
        joined{j} = outerjoin((broken_comp{j}{i}),joined{j},Keys=1,MergeKeys=true);
    end
    joined{j} = table2array(joined{j});
end
toc

%% the data should be ensemble averaged BEFORE being phase averaged.
dg = 1;
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
%%
figure
cmap = jet(8)/1.2;
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
    plot(phase_averaged{i}(:,1),pa_shifted{i},'LineWidth',1,'Color', cmap(i,:))
end
yline(0)
legend('SG-1','SG-2','SG-3','SG-4','SG-5','SG-6','SG-7','SG-8','y = 0')
grid on
xlim([0 360])
xlabel('Degree of rotation')
ylabel('Microstrain')
%title("Phase average at sensor " + i + " with " + degree + " degree ensemble averages")

% figure
% plot(phase_averaged(:,1),phase_averaged(:,2))
% figure
% plot(linspace(-pi,pi,length(phase_averaged(:,1))),phase_averaged(:,2))

% pa_end = phase_averaged(round(1/4*length(phase_averaged)):end,2);
% pa_start = phase_averaged(1:round(1/4*length(phase_averaged)-1),2);
% pa_fixed = [pa_end;pa_start];
% figure
% plot(ens_ave(:,1),pa_fixed)
% grid on
% degree = num2str(dg);
% title("Phase averaged data with " + degree + " degree ensemble averages")

%% averaging only the similar length phases
j=1;
for i = 1:1966
    new_bc{i} = height(broken_comp{3}{i}) > 175 & height(broken_comp{3}{i}) < 300;
end

tic
    njoined = (broken_comp{3}{1});
    for i = 2:height(broken_comp{1})
        if new_bc{i} == 1
            njoined = outerjoin((broken_comp{3}{i}),njoined,Keys=1,MergeKeys=true);
        end
    end
    njoined = table2array(njoined);
toc
%%
dg = 1;
degree = num2str(dg);
nens_ave = [];
for i = 1:183
    for j = 1:length(njoined(:,1))/dg
        nens_ave{i}(j,:) = mean(njoined(dg*(j-1)+1:dg*j,i),'omitnan');
    end
end
nens_ave = cell2mat(nens_ave);
%%
    nphase_averaged = [nens_ave(:,1),mean(nens_ave(:,2:end),2,'omitnan')];

    %subtract the mean
    nmeans = mean(nphase_averaged(:,2),'omitnan');
    npa_a = nphase_averaged(:,2)-nmeans;

    % shifting the zero crossing since its arbitrary
    npa_end = npa_a(round(1/4.5*length(nphase_averaged)):end); 
    npa_start = npa_a(1:round(1/4.5*length(nphase_averaged)-1));
    npa_shifted = [npa_end;npa_start];
    
    hold on
    plot(nphase_averaged(:,1),npa_shifted,'LineWidth',1)

%% histogram of lengths of phases
for i = 1:length(broken_time{1})
    lens(i) = length(broken_time{i});
end
figure
histogram(lens/64,25)
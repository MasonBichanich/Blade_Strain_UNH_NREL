%% This script will read in the ADV data for the tdms files and save it in a .mat structure called DATA.
%gives you the data save in your workspace. The only 'filter' it applies is
%to remove the timestamps (and corresponding data) that are decades off
%Takes about 5 minutes per gig of data depending on machine

clear all

% prompt = 'Which v_link? [1 = 8145, 2 = 28175 ]';
% str = input(prompt,'s');

%path = "..\8.17.23 Deployment\ADV"+str+"\Raw Data\";
path = '../../Data_Files/MODAQ/V_link/8145_Vlink/Raw Data/';
filename = "*.tdms";

% path = '../10.21.22 to 12.10.22 Deployment/Raw Data/VectorVelocityData2';
% filename = '/*_VectorVelocityData2.tdms';

addpath('../tdmsfunctions/')
addpath('../tdmsSubfunctions/')
addpath(path);

unsorted_vlink_data_files = dir(append(path,filename));

table=struct2table(unsorted_vlink_data_files);
sortedtable=sortrows(table,'date');                                 %sorting the files by datenum. the naming is a bit funky
vlink_data_files=table2struct(sortedtable);


for j=1:length(vlink_data_files)                                 

    try                     % this command will ignore the broken files
        DATA(j)=TDMS_readTDMSFile(vlink_data_files(j).name);
    end
   %end
end

%%
N=length(DATA);

Vlink=[];
for j=3:8
    Vlink_dir=[];
    for i=1:N
        if size(DATA(i).data) == [1 9]
            Vlink_cellarray(i,j)=DATA(i).data(j);  
            Vlink_dir=[Vlink_dir; Vlink_cellarray{i,j}'];
        end
        
    end
    Vlink(:,j)=[Vlink_dir];
    j
end
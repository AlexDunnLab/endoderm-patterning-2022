%% CompileAllData %%

% This code takes .txt lists of data (x, y, foxa2, otx2, cdx2) from 
% QuantifyNuclearIntensity.m and compiles them into a sorted list
% containing all data points. For use when performing calculations
% on pooled data, as in Fig. 4A in Cui and Engel et al.

% Written by Kiara Cui
% Last modified December 3, 2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all; rng default;

% prompt for number of files to include
n = inputdlg('How many files?'); n = str2double(n);

allnormdata = []; % initialize matrix to store all data

% extract data from files
for i = 1:n

    [file,path] = uigetfile('*.txt'); % choose '_Area.txt' file
    fname = strcat(path,file); % gets file name from uigetfile
    fname_short = fname(1:end-4);

    fileID = fopen(fname, 'r'); % open file for 'r'eading
    rawdata = readmatrix(fname); % save the numeric data
    fclose(fileID); % close the file after reading

    % normalize data from file by maximum in each dataset (already normalized
    % by Hoescht)
    x = rawdata(:,1);
    normfoxa2 = rawdata(:,3);
    normotx2 = rawdata(:,4); 
    normcdx2 = rawdata(:,5);

    newnormdata = [x, normfoxa2, normotx2, normcdx2];

    % remove outliers from intensity channels
    newnormdata = rmoutliers(newnormdata,'mean'); % nothing from x is removed

    % normalize
    newnormdata(:,2) = newnormdata(:,2)/max(newnormdata(:,2));
    newnormdata(:,3) = newnormdata(:,3)/max(newnormdata(:,3));
    newnormdata(:,4) = newnormdata(:,4)/max(newnormdata(:,4));

    allnormdata = [allnormdata; newnormdata];

end

% sorts by ascending x
allnormdata = sortrows(allnormdata);
data.x = allnormdata(:,1);
data.foxa2 = allnormdata(:,2);
data.otx2 = allnormdata(:,3);
data.cdx2 = allnormdata(:,4);

% choose save destination
folder_name = uigetdir(); cd(folder_name);

%% Save Data

fileID = fopen('EveryDataPoint.txt','w');
writetable(struct2table(data),'EveryDataPoint.txt')
fclose(fileID);
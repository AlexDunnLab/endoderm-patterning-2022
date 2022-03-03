%% AllPointsIntensityPlot %%

% This code takes a .txt file containing all data points (from
% CompileAllData.m) to generate plots with all data points and SD. See 
% Fig. 4A in Cui and Engel et al.

% Written by Kiara Cui
% Last modified January 9, 2022

%% Load Data 

clear; clc; close all; rng default;

[file,path] = uigetfile('*.txt'); % choose '_Area.txt' file
fname = strcat(path,file); % gets file name from uigetfile
fname_short = fname(1:end-4);

fileID = fopen(fname, 'r'); % open file for 'r'eading
rawdata = readmatrix(fname); % save the numeric data
fclose(fileID); % close the file after reading

% label data vectors
x = rawdata(:,1);
foxa2 = rawdata(:,2);
otx2 = rawdata(:,3); 
cdx2 = rawdata(:,4);

%% Plot Preparation 

% choose save destination
folder_name = uigetdir(); cd(folder_name);

% create bins
binedges = (0:25:1000);
[~,binedges,loc] = histcounts(x, binedges); % assign bin indices to each data point
binx = 0.5*(binedges(1:end-1)+binedges(2:end)); % find center of each bin to plot avg value
binx = binx';

%% otx2 plot
binotx2data = accumarray(loc(:), otx2)./accumarray(loc(:),1);
binotx2std = accumarray(loc(:),otx2,[],@std);

% all points scatter plot
figure; 
scatter(x, otx2,'.y');
title('Average OTX2 Intensity vs. x-position')
xlabel('x (\mum)'); ylabel('Average OTX2 intensity')
xlim([0 1000]); ylim([0 1])
saveas(gcf,'OTX2 Scatter','epsc'); saveas(gcf,'OTX2 Scatter','tiffn');

% average and SD
figure; 
errorbar(binx, binotx2data, binotx2std,'-k');
title('Average OTX2 Intensity vs. x-position')
xlabel('x (\mum)'); ylabel('Average OTX2 intensity')
xlim([0 1000]); ylim([0 1])
saveas(gcf,'OTX2 SD','epsc'); saveas(gcf,'OTX2 SD','tiffn');

%% cdx2 plot
bincdx2data = accumarray(loc(:), cdx2)./accumarray(loc(:),1);
bincdx2std = accumarray(loc(:),cdx2,[],@std);

% all points scatter plot
figure; 
scatter(x, cdx2,'.m');
title('Average CDX2 Intensity vs. x-position')
xlabel('x (\mum)'); ylabel('Average CDX2 intensity')
xlim([0 1000]); ylim([0 1])
% save in vector and tif format
saveas(gcf,'CDX2 Scatter','epsc'); saveas(gcf,'CDX2 Scatter','tiffn');

% average and SD
figure; 
errorbar(binx, bincdx2data, bincdx2std,'-k');
title('Average CDX2 Intensity vs. x-position')
xlabel('x (\mum)'); ylabel('Average CDX2 intensity')
xlim([0 1000]); ylim([0 1])
% save in vector and tif format
saveas(gcf,'CDX2 SD','epsc'); saveas(gcf,'CDX2 SD','tiffn');
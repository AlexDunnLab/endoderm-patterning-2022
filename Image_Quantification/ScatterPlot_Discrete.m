%% Generating Pooled Mock FACS Plots %%

% This code takes lists of data (x, y, foxa2, otx2, cdx2) from discrete
% experiments and compiles all data points into a pooled scatter plot 
% comparing any two intensities, with the goal of identifying distinct populations.
% See Figure 1H in Cui and Engel et al.

% Written by Kiara Cui
% Last modified September 20, 2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all; rng default;

% prompt for number of files to compile
n = inputdlg('How many files?'); n = str2double(n);

allnormdata = []; % initialize matrix to store all data

% no need for additional normalization
maxotx2 = 1;
maxcdx2 = 1; 
maxfoxa2 = 1;

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
normfoxa2 = rawdata(:,3)/maxfoxa2; % normalize by own foxa2
normotx2 = rawdata(:,4)/maxotx2; 
normcdx2 = rawdata(:,5)/maxcdx2;

newnormdata = [normfoxa2, normotx2, normcdx2];

% remove outliers from any channel
newnormdata = rmoutliers(newnormdata);

allnormdata = [allnormdata; newnormdata];

end

%% Generate Plots

% choose save destination
folder_name = uigetdir(); cd(folder_name);

% foxa2 histogram
figure; histogram(allnormdata(:,1));
xlabel('FOXA2'); ylabel('Count');
saveas(gcf,'Foxa2 Histogram','tiffn');
saveas(gcf,'Foxa2 Histogram','epsc'); 

% otx2 histogram
figure; histogram(allnormdata(:,2));
xlabel('OTX2'); ylabel('Count');
saveas(gcf,'Otx2 Histogram','tiffn');
saveas(gcf,'Otx2 Histogram','epsc'); 

% cdx2 histogram
figure; histogram(allnormdata(:,3));
xlabel('CDX2'); ylabel('Count');
saveas(gcf,'Cdx2 Histogram','tiffn');
saveas(gcf,'Cdx2 Histogram','epsc'); 

% Scatter plot for otx2 vs cdx2

figure; plot(allnormdata(:,2), allnormdata(:,3),'.','MarkerSize',1);
xlabel('OTX2'); ylabel('CDX2');
xlim([0,2]); ylim([0,2]);
% save
saveas(gcf,'Mock FACS for Otx2 vs Cdx2','tiffn');
saveas(gcf,'Mock FACS for Otx2 vs Cdx2','epsc'); 

% Scatter plot for foxa2 vs otx2

figure; plot(allnormdata(:,1), allnormdata(:,2),'.','MarkerSize',1)
xlabel('FOXA2'); ylabel('OTX2');
xlim([0,2]); ylim([0,2]);
% save
saveas(gcf,'Mock FACS for Foxa2 vs Otx2','tiffn');
saveas(gcf,'Mock FACS for Foxa2 vs Otx2','epsc'); 

% Scatter plot for foxa2 vs cdx2

figure; plot(allnormdata(:,1), allnormdata(:,3),'.','MarkerSize',1)
xlabel('FOXA2'); ylabel('CDX2');
xlim([0,2]); ylim([0,2]);
% save in tif format
saveas(gcf,'Mock FACS for Foxa2 vs Cdx2','tiffn');
saveas(gcf,'Mock FACS for Foxa2 vs Cdx2','epsc'); 

%% Save Data

pooleddata.allnormdata = allnormdata;
fileID = fopen('PooledRawData.txt','w');
writetable(struct2table(pooleddata),'PooledRawData.txt')
fclose(fileID);

function y = ceiln(x,n)
    y = ceil(x*10^n)/(10^n);
end

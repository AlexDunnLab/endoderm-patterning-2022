%% Generating Pooled Scatter Plots %%

% This code takes lists of raw data (x, y, foxa2, otx2, cdx2) from gradient
% experiments and compiles all data points into a pooled scatter plot 
% comparing any two intensities, with the goal of identifying distinct populations.
% See Figure 1G in Cui and Engel et al.

% Written by Kiara Cui
% Last modified August 20, 2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all; rng default;

% prompt for number of files to incorporate
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
normfoxa2 = rawdata(:,3); %/max(rawdata(:,3));
normotx2 = rawdata(:,4); %/max(rawdata(:,4)); 
normcdx2 = rawdata(:,5);% /max(rawdata(:,5));

newnormdata = [normfoxa2, normotx2, normcdx2];

% remove outliers from any channel

newnormdata = rmoutliers(newnormdata,'mean');

newnormdata(:,1) = newnormdata(:,1)/max(newnormdata(:,1));
newnormdata(:,2) = newnormdata(:,2)/max(newnormdata(:,2));
newnormdata(:,3) = newnormdata(:,3)/max(newnormdata(:,3));

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

% Scatter plot for otx2 vs cdx2 intensities

figure; plot(allnormdata(:,2), allnormdata(:,3),'.','MarkerSize',1);
xlabel('OTX2'); ylabel('CDX2');
if max(allnormdata(:,2)) > max(allnormdata(:,3))
    limit1 = ceiln(max(allnormdata(:,2)),1); % set limit to one decimal place above maximum
    xlim([0, limit1]); ylim([0, limit1]);
elseif max(allnormdata(:,3)) > max(allnormdata(:,2))
    limit1 = ceiln(max(allnormdata(:,3)),1); % set limit to one decimal place above maximum
    xlim([0, limit1]); ylim([0, limit1]);
end
% save
saveas(gcf,'Mock FACS for Otx2 vs Cdx2','tiffn');
saveas(gcf,'Mock FACS for Otx2 vs Cdx2','epsc'); 

% Scatter plo for foxa2 vs otx2

figure; plot(allnormdata(:,1), allnormdata(:,2),'.','MarkerSize',1)
xlabel('FOXA2'); ylabel('OTX2');
if max(allnormdata(:,1)) > max(allnormdata(:,2))
    limit2 = ceiln(max(allnormdata(:,1)),1); % set limit to one decimal place above maximum
    xlim([0, limit2]); ylim([0, limit2]);
elseif max(allnormdata(:,2)) > max(allnormdata(:,1))
    limit2 = ceiln(max(allnormdata(:,2)),1); % set limit to one decimal place above maximum
    xlim([0, limit2]); ylim([0, limit2]);
end
% save
saveas(gcf,'Mock FACS for Foxa2 vs Otx2','tiffn');
saveas(gcf,'Mock FACS for Foxa2 vs Otx2','epsc'); 

% Scatter plot for foxa2 vs cdx2

figure; plot(allnormdata(:,1), allnormdata(:,3),'.','MarkerSize',1)
xlabel('FOXA2'); ylabel('CDX2');
if max(allnormdata(:,1)) > max(allnormdata(:,3))
    limit3 = ceiln(max(allnormdata(:,1)),1); % set limit to one decimal place above maximum
    xlim([0, limit3]); ylim([0, limit3]);
elseif max(allnormdata(:,3)) > max(allnormdata(:,1))
    limit3 = ceiln(max(allnormdata(:,3)),1); % set limit to one decimal place above maximum
    xlim([0, limit3]); ylim([0, limit3]);
end
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

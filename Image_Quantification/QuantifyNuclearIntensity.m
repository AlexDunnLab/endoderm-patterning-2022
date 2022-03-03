%% QuantifyNuclearIntensity.m %%

% Adapted from framework from Christina Hueschen in 04/2020 
% Last modified by Kiara Cui on 3/15/21

% This script analyzes multichannel .tif confocal images by quantifying
% and plotting intensites from four channels: 
% (1) Hoechst (405 nm, cyan), 
% (2) FOXA2 (488 nm, green), 
% (3) CDX2 (555 nm, magenta), and
% (4) OTX2 (647 nm, yellow).
% Output .txt files contain data that are further processed to generate
% final graphs in Cui and Engel et al.

%% Set parameters and load images

clearvars % clear all variables
clc % clear matlab workspace
close all % close all existing windows
warning('off','all'); % Suppress warnings for faster exec.

% Set parameters
pixSize = 1.6605; % microns per pixel, defined in image metadata

% Load the images: select folder containing images
folder_name = uigetdir(); cd(folder_name);

% this struct saves info about the images to be analyzed
dirFileList = dir([folder_name, '/SUM*.tif']); % only reads in TIF files whose names start with SUM

%% Select areas for analysis
 
im = imread(dirFileList(1).name, 1); % load channel 1 (Hoechst) and save in as matrix im--will be used for segmentation
[imsize_y, imsize_x] = size(im); % find image size

% Apply a median filter to reduce the noise in the image. 
imFilt = medfilt2(im); % sets each pixel to median pixel of the 3x3 box centered there
imNorm = mat2gray(imFilt); % scales values (intensity) so that min = 0 and max = 1
    
% Choose left and right boundaries for channel / cell area
trial = 0; % initialize
while trial >= 0
    % Find edges of cell sheet
    figure; imshow(im, []);
    message = sprintf('After pressing okay: Click to define the left edge of the cell sheet, then the right edge');
    uiwait(msgbox(message));
    [cellsX, cellsY] = ginput(2); % Use matlab function ginput to save x and y coordinates of your click
   
    % save coordinates as variables
    leftedgecells = cellsX(1);
    rightedgecells = cellsX(2);
    
    % display selected area so user can confirm
    hold on
    line([leftedgecells leftedgecells],[0 imsize_y],'Color','cyan','LineStyle','--')
    line([rightedgecells rightedgecells],[0 imsize_y],'Color','cyan','LineStyle','--')
    hold off
    
    % confirm selected area
    areacheck = questdlg('Corresponds to the marked area on the image?',...
        'Area check','Yes','No','Yes');
    switch areacheck
        case 'Yes'
            break % Quit while-loop
        case 'No'
            trial = trial + 1; % Redo while-loop
    end
end

% select left and right edge of cell channel
while trial >= 0
   
    message = sprintf('After pressing okay, click to define the left edge of the channel. The right edge will be calculated then shown as well.');
    uiwait(msgbox(message));
    % Use matlab function ginput to save x and y coordinates of your click
    [channelX, channelY] = ginput(1);
    
    % save channel edge coordinates; automatically place right channel edge
    leftchanneledge = channelX(1);
    rightchanneledge = channelX(1) + 1000/pixSize;
    
    % display selected area so user can confirm
    hold on
    line([leftchanneledge leftchanneledge],[0 imsize_y],'Color','red','LineStyle','--')
    line([rightchanneledge rightchanneledge],[0 imsize_y],'Color','red','LineStyle','--')
    hold off
    
    % confirm selected area
    areacheck = questdlg('Corresponds to the marked area on the image?',...
        'Area check','Yes','No','Yes');
    switch areacheck
        case 'Yes'
            break % Quit while-loop
        case 'No'
            trial = trial + 1; % Redo while-loop
    end
end

% select top and bottom of analysis area, avoid areas with incomplete cell
% coverage and edge effects
while trial >= 0
   
    message = sprintf('After pressing okay, click to define the top and bottom of the area to be analyzed.');
    uiwait(msgbox(message));
    % Use matlab function ginput to save x and y coordinates of your click
    [bandX, bandY] = ginput(2);
    
    % save band edge coordinates
    bandtop = bandY(1);
    bandbottom = bandY(2);
    
    % display selected area(s) so user can confirm
    hold on
    line([0 imsize_x],[bandtop bandtop],'Color','green','LineStyle','--')
    line([0 imsize_x],[bandbottom bandbottom],'Color','green','LineStyle','--')
    hold off
    
    % confirm selected area
    areacheck = questdlg('Corresponds to the marked area on the image?',...
        'Area check','Yes','No','Yes');
    switch areacheck
        case 'Yes'
            break % Quit while-loop
        case 'No'
            trial = trial + 1; % Redo while-loop
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% Choose a region to use to take a background intensity measurement
message = sprintf('Define an ROI for a background measurement (ideally within the channel but where there are no cells). After pressing Okay, first click to choose the ROI vertices.\nClose the polygon by clicking on your first point.\nDouble-click inside polygon to accept it.');
uiwait(msgbox(message));
% Use roipoly function to create mask from user chosen polygon
BGmask = roipoly(imNorm); 
    
% Try a threshold
thresh = adaptthresh(imNorm, 0.75); 
% "The graythresh function uses Otsu's method, which chooses the threshold
% to minimize the intraclass variance of the black and white pixels"
% second parameter is sensitivity, default 0.5, higher means include
% more as foreground
imThresh = imNorm > thresh; %keep only pixels above threshold
    
%% Segment nuclei and save as objects

% quasi-euclidean distance transform image; assigns values based on nearest 
% nonzero neighbor such that only nucleus centers are preserved
D = -bwdist(~imThresh,'quasi-euclidean'); figure; imshow(D,[]) % show resulting image

% compute watershed transform and use ridge lines to segment nuclei
DL = watershed(D); imThresh2 = imThresh; imThresh2(DL == 0) = 0; figure; imshow(imThresh2) % show resulting image

% filter out tiny local minima - should generate tiny dot in the center of
% each nucleus
mask = imextendedmin(D,0.5,4); figure; imshowpair(imThresh,mask,'blend') % show resulting image

% modify distance transform so minima only occur at desired locations
D2 = imimposemin(D, mask); LD2 = watershed(D2); imThresh3 = imThresh;
imThresh3(LD2 == 0) = 0; figure; imshow(imThresh3) % show resulting image

% identify nuclei edges by where gradient of intensity is a maximum
gmag = imgradient(imThresh3); figure; imshow(gmag,[]) % show edges in figure

% implement watershed segmentation based on gradient image
imWS = watershed(gmag,4); Lrgb = label2rgb(imWS,'spring','b','shuffle'); % convert to RGB image for visualization %     figure; imshow(Lrgb) % display segmentation

% Label non-zero pixel clusters as objects
imLab = bwlabel(imWS,4); % label connected components in 2D binary image

% Now we want to throw out all objects that don't meet our criteria: a nucleus-sized area and shape, 
% and a location within the boundaries of our channel 
clear props
% use regionprops to find area, solidity (convex area), centroid, circularity of a binary array, imLab
props = regionprops(imLab, 'Area', 'Solidity', 'Centroid','Circularity');
approvedObj = zeros(size(im));

for i = 1:length(props)
    if props(i).Area > 15 && props(i).Area < 1250 && props(i).Solidity > 0.7 ...
            && props(i).Centroid(1) > leftedgecells && props(i).Centroid(1) < rightedgecells...
            && props(i).Centroid(2) < bandbottom && props(i).Centroid(2) > bandtop...
            && props(i).Circularity > 0.7
        % keeps only objects of certain size, circularity, and position from imLab
        approvedObj = approvedObj + (imLab==i); 
    end
end

% relabel only these approved objects
imLab = bwlabel(approvedObj); % used to be bwlabel
boundaries = bwboundaries(approvedObj,'noholes'); % save boundaries of each cell

%% Measure intensities in each channel (Hoechst, FOXA2, CDX2, OTX2)

% Save areas and positions of segmented nuclei
clear props
% use regionprops to find area and centroid position
props = regionprops(imLab, 'Area', 'Centroid');

% Measure average Hoechst intensities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear im
im = imread(dirFileList(1).name,1); %pull in Hoechst image
BGimHoechst = double(im).*BGmask;
BGimHoechst(BGimHoechst == 0) = NaN;
BGHoechst = nanmean(nanmean(BGimHoechst,1),2);

segHoechst = double(im).*approvedObj;
segHoechst(segHoechst == 0) = NaN;

clear Hoechstprops
clear HoechstrawInt; HoechstrawInt = [];
clear HoechstintMinusBG; HoechstintMinusBG = [];

Hoechstprops = regionprops(imLab, im, 'MeanIntensity');
tempRawInt = cell2mat({Hoechstprops(1:end).MeanIntensity});
tempIntMinusBG = tempRawInt - BGHoechst;
HoechstrawInt = [HoechstrawInt tempRawInt];
HoechstintMinusBG = [HoechstintMinusBG tempIntMinusBG];

% Measure average FOXA2 intensities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

im = imread(dirFileList(1).name,2); % pull in foxa2 image 
BGimFoxa2 = double(im).*BGmask;
BGimFoxa2(BGimFoxa2 == 0) = NaN;
BGFoxa2 = nanmean(nanmean(BGimFoxa2,1),2); % finds mean excluding NaN values
segFoxa2 = double(im).*approvedObj;
segFoxa2(segFoxa2 == 0) = NaN;

clear Foxa2props
clear Foxa2rawInt; Foxa2rawInt = [];
clear Foxa2intMinusBG; Foxa2intMinusBG = [];

Foxa2props = regionprops(imLab, im, 'MeanIntensity');
tempRawInt = cell2mat({Foxa2props(1:end).MeanIntensity});
tempIntMinusBG = tempRawInt - BGFoxa2;
Foxa2rawInt = [Foxa2rawInt tempRawInt];
Foxa2intMinusBG = [Foxa2intMinusBG tempIntMinusBG];

% Measure CDX2 intensities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear im
im = imread(dirFileList(1).name,3); %pull in cdx2 image
BGimCdx2 = double(im).*BGmask;
BGimCdx2(BGimCdx2 == 0) = NaN;
BGCdx2 = nanmean(nanmean(BGimCdx2,1),2);

segCdx2 = double(im).*approvedObj;
segCdx2(segCdx2 == 0) = NaN;

clear Cdx2props
clear Cdx2rawInt; Cdx2rawInt = [];
clear Cdx2intMinusBG; Cdx2intMinusBG = [];

Cdx2props = regionprops(imLab, im, 'MeanIntensity');
tempRawInt = cell2mat({Cdx2props(1:end).MeanIntensity});
tempIntMinusBG = tempRawInt - BGCdx2;
Cdx2rawInt = [Cdx2rawInt tempRawInt];
Cdx2intMinusBG = [Cdx2intMinusBG tempIntMinusBG];

% Measure OTX2 intensities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear im
im = imread(dirFileList(1).name,4); %pull in Otx2 image
BGimOtx2 = double(im).*BGmask;
BGimOtx2(BGimOtx2 == 0) = NaN;
BGOtx2 = nanmean(nanmean(BGimOtx2,1),2);

segOtx2 = double(im).*approvedObj;
segOtx2(segOtx2 == 0) = NaN;

clear Otx2props
clear Otx2rawInt; Otx2rawInt = [];
clear Otx2intMinusBG; Otx2intMinusBG = [];

Otx2props = regionprops(imLab, im, 'MeanIntensity');
tempRawInt = cell2mat({Otx2props(1:end).MeanIntensity});
tempIntMinusBG = tempRawInt - BGOtx2;
Otx2rawInt = [Otx2rawInt tempRawInt];
Otx2intMinusBG = [Otx2intMinusBG tempIntMinusBG];

% Save data into matrix %%%%%%%%%%%%%%%%%%
data = zeros(size(props,1),1);

for i = 1:size(props,1) 
    % save areas in first column of array
    data(i,1) = props(i).Area * pixSize^2; 
    % save x positions. Set edge of channel as x=0! Convert to microns
    data(i,2) = (props(i).Centroid(1) - leftchanneledge) * pixSize;
    % save y position in microns
    data(i,3) = props(i).Centroid(2) * pixSize;
end
    
% save intensities
data(:,4) = Foxa2intMinusBG';
data(:,5) = Cdx2intMinusBG';
data(:,6) = Otx2intMinusBG';
data(:,7) = HoechstintMinusBG'; 

close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Discard outliers based on lack of nuclear signal

medHoechst_unfiltered = median(data(:,7)); % find current median
norm_Hoechst = data(:,7)./medHoechst_unfiltered; % use Hoechst normalized by median to identify segmented objects that have no nuclear signal
outliers = find(norm_Hoechst < 0.3); % find indices of outliers (normalized Hoechst intensity < 0.3)
for i = size(outliers):-1:1 % loop from bottom up since size changes
    data(outliers(i),:) = []; % remove from data matrix
    boundaries(outliers(i),:) = []; % remove from boundaries matrix
end

% Simplify name of position data
areadata = data(:,1);
xdata = data(:,2);
ydata = data(:,3);
avgFoxa2 = data(:,4);
avgCdx2 = data(:,5);
avgOtx2 = data(:,6);
avgHoechst = data(:,7);

% Correct negative intensities (due to background subtraction)
for n = 1:length(xdata)
    if avgHoechst(n) < 0
        avgHoechst(n) = 0;
    end
    if avgFoxa2(n) < 0
        avgFoxa2(n) = 0;
    end
    if avgCdx2(n) < 0
        avgCdx2(n) = 0;
    end
    if avgOtx2(n) < 0
        avgOtx2(n) = 0;
    end
end

%% Create bins
  
% create structure to store output
output = struct('x',zeros(40,1),'count',zeros(40,1),'density',zeros(40,1),...
    'Hoechstdata',zeros(40,1),'Hoechststdev',zeros(40,1),...
    'foxa2data',zeros(40,1),'foxa2stdev',zeros(40,1),...
    'cdx2data',zeros(40,1),'cdx2stdev',zeros(40,1),...
    'otx2data',zeros(40,1),'otx2stdev',zeros(40,1));

binedges = (-500:50:1500);
[~,binedges,loc] = histcounts(xdata, binedges); % assign bin indices to each data point
output.x = 0.5*(binedges(1:end-1)+binedges(2:end)); % find center of each bin to plot avg value
output.x = output.x';
    
%% cell density vs. position

% define bins and calculate area of each one
binheight = (bandbottom - bandtop)*pixSize; binwidth = 50; % microns
output.count = accumarray(loc(:),1);
output.density = output.count*(10^6)/(binheight*binwidth);

% backfill entries where no cells are present
for k = length(output.count)+1:40
    output.density(k) = 0;
    output.count(k) = 0;
end
% plot 
figure; plot(output.x,output.density,'-','LineWidth',1)
title('Nuclear density'); xlabel('x (\mum)'); ylabel('Nuclei per mm^2')
xlim([-500 1500]); ylim([0 3500]) 
% save in vector format
saveas(gcf,'Nuclear Density','epsc')
% save in tif format
saveas(gcf,'Nuclear Density','tiffn')
    
%% Hoechst intensity vs. position
    
medavgHoechst = median(avgHoechst);
output.Hoechstdata = accumarray(loc(:), avgHoechst./medavgHoechst)./accumarray(loc(:),1);
output.Hoechststdev = accumarray(loc(:),avgHoechst./medavgHoechst,[],@std);

for k = length(output.Hoechstdata)+1:1:40
    output.Hoechstdata(k) = NaN;
    output.Hoechststdev(k) = NaN;
end

figure; hold on
scatter(xdata, avgHoechst./medavgHoechst,'.');
shadedErrorBar(output.x, output.Hoechstdata, output.Hoechststdev,'lineProps','-b','transparent',1);
hold off
title('Average Hoechst Intensity vs. x-position')
xlabel('x (\mum)'); ylabel('Average Hoechst intensity, normalized by median')
% legend('Raw Data','Binned Averages','Smoothed Averages')
xlim([-500 1500]); ylim([0 inf])
% save in vector and tif format
saveas(gcf,'Hoechst Intensity','epsc'); saveas(gcf,'Hoechst Intensity','tiffn');

% make separate plot scaled by maximum
figure
shadedErrorBar(output.x, output.Hoechstdata/max(output.Hoechstdata),...
    output.Hoechststdev/max(output.Hoechstdata),'lineProps','-b','transparent',1);
title('Normalized Hoechst Intensity')
xlabel('x (\mum)'); ylabel('Normalized Hoechst Intensity'); xlim([-500 1500]); ylim([0 1.5])
% save in vector and tif format
saveas(gcf,'Normalized Hoechst Intensity','epsc'); saveas(gcf,'Normalized Hoechst Intensity','tiffn');
    
%% foxa2 intensity vs. position

output.foxa2data = accumarray(loc(:), avgFoxa2./avgHoechst)./accumarray(loc(:),1);
output.foxa2stdev = accumarray(loc(:),avgFoxa2./avgHoechst,[],@std);

for k = length(output.foxa2data)+1:1:40
    output.foxa2data(k) = NaN;
    output.foxa2stdev(k) = NaN;
end

figure; hold on
scatter(xdata, avgFoxa2./avgHoechst,'.');
shadedErrorBar(output.x, output.foxa2data, output.foxa2stdev,'lineProps','-g','transparent',1);
hold off
title('Average FOXA2 Intensity vs. x-position')
xlabel('x (\mum)'); ylabel('Average FOXA2 intensity, normalized by Hoechst')
% legend('Raw Data','Binned Averages','Smoothed Averages')
xlim([-500 1500]); ylim([0 inf])
% save in vector and tif format
saveas(gcf,'FOXA2 Intensity','epsc'); saveas(gcf,'FOXA2 Intensity','tiffn');

% make separate plot scaled by maximum
figure
shadedErrorBar(output.x, output.foxa2data/max(output.foxa2data),...
    output.foxa2stdev/max(output.foxa2data),'lineProps','-g','transparent',1);
title('Normalized FOXA2 Intensity')
xlabel('x (\mum)'); ylabel('Normalized FOXA2 Intensity'); xlim([-500 1500]); ylim([0 1.5])
% save in vector and tif format
saveas(gcf,'Normalized FOXA2 Intensity','epsc'); saveas(gcf,'Normalized FOXA2 Intensity','tiffn');

%% otx2 intensity vs. position
 
output.otx2data = accumarray(loc(:), avgOtx2./avgHoechst)./accumarray(loc(:),1);
output.otx2stdev = accumarray(loc(:),avgOtx2./avgHoechst,[],@std);

for k = length(output.otx2data)+1:1:40
    output.otx2data(k) = NaN;
    output.otx2stdev(k) = NaN;
end

figure; hold on
scatter(xdata, avgOtx2./avgHoechst,'.');
shadedErrorBar(output.x, output.otx2data, output.otx2stdev,'lineProps','-y','transparent',1);
hold off
title('Average OTX2 Intensity vs. x-position')
xlabel('x (\mum)'); ylabel('Average OTX2 intensity, normalized by Hoechst')
% legend('Raw Data','Binned Averages','Smoothed Averages')
xlim([-500 1500]); ylim([0 inf])
% save in vector and tif format
saveas(gcf,'OTX2 Intensity','epsc'); saveas(gcf,'OTX2 Intensity','tiffn');

% make separate plot scaled by maximum
figure
shadedErrorBar(output.x, output.otx2data/max(output.otx2data),...
    output.otx2stdev/max(output.otx2data),'lineProps','-y','transparent',1);
title('Normalized OTX2 Intensity')
xlabel('x (\mum)'); ylabel('Normalized OTX2 Intensity'); xlim([-500 1500]); ylim([0 1.5])
% save in vector and tif format
saveas(gcf,'Normalized OTX2 Intensity','epsc'); saveas(gcf,'Normalized OTX2 Intensity','tiffn');

%% cdx2 intensity vs. position

output.cdx2data = accumarray(loc(:), avgCdx2./avgHoechst)./accumarray(loc(:),1);
output.cdx2stdev = accumarray(loc(:),avgCdx2./avgHoechst,[],@std);

for k = length(output.cdx2data)+1:1:40
    output.cdx2data(k) = NaN;
    output.cdx2stdev(k) = NaN;
end

figure; hold on
scatter(xdata, avgCdx2./avgHoechst,'.');
shadedErrorBar(output.x, output.cdx2data, output.cdx2stdev,'lineProps','-m','transparent',1);
hold off
title('Average CDX2 Intensity vs. x-position')
xlabel('x (\mum)'); ylabel('Average CDX2 intensity, normalized by Hoechst')
% legend('Raw Data','Binned Averages','Smoothed Averages')
xlim([-500 1500]); ylim([0 inf])
% save in vector and tif format
saveas(gcf,'CDX2 Intensity','epsc'); saveas(gcf,'CDX2 Intensity','tiffn');

% make separate plot scaled by maximum
figure
shadedErrorBar(output.x, output.cdx2data/max(output.cdx2data),...
    output.cdx2stdev/max(output.cdx2data),'lineProps','-m','transparent',1);
title('Normalized CDX2 Intensity')
xlabel('x (\mum)'); ylabel('Normalized CDX2 Intensity'); xlim([-500 1500]); ylim([0 1.5])
% save in vector and tif format
saveas(gcf,'Normalized CDX2 Intensity','epsc'); saveas(gcf,'Normalized CDX2 Intensity','tiffn');

%% double positive?
% check for cells that are both OTX2 and CDX2 positive, defined by a
% user-specified threshold

normCdx2 = (avgCdx2./avgHoechst)/max(output.cdx2data);
normOtx2 = (avgOtx2./avgHoechst)/max(output.otx2data);
Cdx2Otx2Ratio = normCdx2./normOtx2;

figure; scatter(normCdx2, normOtx2)
xlabel('CDX2'); ylabel('OTX2');
xlim([0, 1.2]); ylim([0, 1.2]);

% mark double positive cells
similarratio = find(Cdx2Otx2Ratio > 0.75 & Cdx2Otx2Ratio < 1.25);
realpositiveCdx2 = find(normCdx2 > 0.4);
realpositiveOtx2 = find(normOtx2 > 0.4);
realpositives = intersect(realpositiveCdx2, realpositiveOtx2);
doublepositives = intersect(realpositives, similarratio);

figure; imfoxa2 = imread(dirFileList(1).name,2); imshow(imfoxa2,[])
hold on
for k = 1:length(realpositiveCdx2)
    boundary = boundaries{realpositiveCdx2(k)};
    plot(boundary(:,2), boundary(:,1),'m')
end
for k = 1:length(realpositiveOtx2)
    boundary = boundaries{realpositiveOtx2(k)};
    plot(boundary(:,2), boundary(:,1),'y')
end

if isempty(doublepositives) == 0
    for k = 1:length(doublepositives)
        boundary = boundaries{doublepositives(k)};
        plot(boundary(:,2), boundary(:,1),'w')
    end
    hold off
elseif isempty(doublepositives) == 1
    disp('There are no cells that were positive for OTX2 and CDX2.')
end

saveas(gcf,'Boundaries','tiffn');
    
%% Write data to text file

fileID = fopen('Analysis.txt','w');
writetable(struct2table(output),'Analysis.txt')
fclose(fileID);

%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout=shadedErrorBar(x,y,errBar,varargin)
% Rob Campbell - November 2009

% Parse input arguments
narginchk(3,inf)

params = inputParser;
params.CaseSensitive = false;
params.addParameter('lineProps', '-k', @(x) ischar(x) | iscell(x));
if (sum( size(ver('MATLAB'))) > 0  )
  params.addParameter('transparent', true, @(x) islogical(x) || x==0 || x==1);
elseif (sum( size(ver('Octave'))) > 0  )
  params.addParameter('transparent', false, @(x) islogical(x) || x==0 || x==1);
end
params.addParameter('patchSaturation', 0.2, @(x) isnumeric(x) && x>=0 && x<=1);

params.parse(varargin{:});

%Extract values from the inputParser
lineProps =  params.Results.lineProps;
transparent =  params.Results.transparent;
patchSaturation = params.Results.patchSaturation;

if ~iscell(lineProps), lineProps={lineProps}; end

%Process y using function handles if needed to make the error bar dynamically
if iscell(errBar) 
    fun1=errBar{1}; fun2=errBar{2}; errBar=fun2(y); y=fun1(y);
else
    y=y(:).';
end

if isempty(x)
    x=1:length(y);
elseif sum( size(ver('MATLAB'))) > 0 
    x=x(:).';
end

%Make upper and lower error bars if only one was specified
if length(errBar)==length(errBar(:))
    errBar=repmat(errBar(:)',2,1);
else
    s=size(errBar); f=find(s==2);
    if isempty(f), error('errBar has the wrong size'), end
    if f==2, errBar=errBar'; end
end

% Check for correct x, errbar formats
x_size = size(x);

if (length(x) ~= length(errBar) && sum( size(ver('MATLAB'))) > 0 )
    error('length(x) must equal length(errBar)')
elseif( ( length(x) ~= length(errBar) && checkOctave_datestr(x) == false ) ...
            && sum( size(ver('Octave'))) > 0  )
    error('length(x) must equal length(errBar) or x must have valid datestr')
end

% Log the hold status so we don't change
initialHoldStatus=ishold;
if ~initialHoldStatus, hold on,  end

H = makePlot(x,y,errBar,lineProps,transparent,patchSaturation);

if ~initialHoldStatus, hold off, end

if nargout==1
    varargout{1}=H;
end
end

function H = makePlot(x,y,errBar,lineProps,transparent,patchSaturation)
    % Determine host application
    if (sum( size(ver('MATLAB'))) > 0  )
      hostName = 'MATLAB';
    elseif (sum(size(ver('Octave'))) > 0)
      hostName = 'Octave';
    end % if  
    % Plot to get the parameters of the line
    if hostName == 'MATLAB'
      H.mainLine=plot(x,y,lineProps{:}); 
    elseif hostName == 'Octave'
      boolxDatestr = checkOctave_datestr(x);
      if boolxDatestr
        x = datenum(x);
        x = x(:).';
        H.mainLine=plot(x,y,lineProps{:});
        datetick(gca);
      else
        H.mainLine=plot(x,y,lineProps{:});
      end
    end
    % Tag the line so we can easily access it
    H.mainLine.Tag = 'shadedErrorBar_mainLine';
    % Work out the color of the shaded region and associated lines.
    % Here we have the option of choosing alpha or a de-saturated
    % solid colour for the patch surface.
    mainLineColor=get(H.mainLine,'color');
    edgeColor=mainLineColor+(1-mainLineColor)*0.55;
    if transparent
        faceAlpha=patchSaturation;
        patchColor=mainLineColor;
    else
        faceAlpha=1;
        patchColor=mainLineColor+(1-mainLineColor)*(1-patchSaturation);
    end
    % Calculate the error bars
    uE=y+errBar(1,:);lE=y-errBar(2,:);
    % Make the patch (the shaded error bar)
    yP=[lE,fliplr(uE)]; xP=[x,fliplr(x)];
    %remove nans otherwise patch won't work
    xP(isnan(yP))=[]; yP(isnan(yP))=[];
    if isdatetime(x) && strcmp(hostName,'MATLAB')
      H.patch=patch(datenum(xP),yP,1);
    else
      H.patch=patch(xP,yP,1);
    end
    set(H.patch,'facecolor',patchColor, ...
        'edgecolor','none', ...
        'facealpha',faceAlpha, ...
        'HandleVisibility', 'off', ...
        'Tag', 'shadedErrorBar_patch')
    %Make pretty edges around the patch. 
    H.edge(1)=plot(x,lE,'-'); H.edge(2)=plot(x,uE,'-');
    set([H.edge], 'color',edgeColor, 'HandleVisibility','off','Tag', 'shadedErrorBar_edge')
    % Ensure the main line of the plot is above the other plot elements
    if hostName == 'MATLAB'
      if strcmp(get(gca,'YAxisLocation'),'left') %Because re-ordering plot elements with yy plot is a disaster
        uistack(H.mainLine,'top')
      end
    elseif hostName == 'Octave'
      % create the struct from scratch by temp.
      H = struct('mainLine', H.mainLine, ...
      'patch', H.patch, ...
      'edge', H.edge);
    end
end

%% Simple try/catch for casting datenums, requireing valid datestr
function boolDate = checkOctave_datestr(x)
  boolDate = true;
  try
    datenum(x)
  catch
    boolDate = false;
  end
end
%%%%%%%%%%%%

%% IF Data Compiler for Project Psyche Data %%

% Generates average data and plots for analyzed immunofluorescence data
% (Fig. 1F) in Cui and Engel et al.

% Written by Kiara Cui
% Last updated April 5, 2021

close all;

% prompt for number of files to average
n = inputdlg('How many files to average?'); n = str2double(n);
x = linspace(-475,1475,40); x = x';

% Initialize matrices for averaging
count = zeros(40, n);
density = zeros(40, n);
hoeschtdata = zeros(40, n);
hoeschtstdev = zeros(40, n);
foxa2data = zeros(40, n);
foxa2stdev = zeros(40, n);
cdx2data = zeros(40, n);
cdx2stdev = zeros(40, n);
otx2data = zeros(40, n);
otx2stdev = zeros(40, n);

% extract data from files
for i = 1:n

[file,path] = uigetfile('*.txt'); % choose '_Area.txt' file
fname = strcat(path,file); % gets file name from uigetfile
fname_short = fname(1:end-4);

fileID = fopen(fname, 'r'); % open file for 'r'eading
alldata = readmatrix(fname); % save the numeric data
fclose(fileID); % close the file after reading

% grab columns from file
count(:,i) = alldata(:,2);
density(:,i) = alldata(:,3);
hoeschtdata(:,i) = alldata(:,4);
hoeschtstdev(:,i) = alldata(:,5);
foxa2data(:,i) = alldata(:,6);
foxa2stdev(:,i) = alldata(:,7);
cdx2data(:,i) = alldata(:,8);
cdx2stdev(:,i) = alldata(:,9);

% ask whether to include otx2 data (some files include cdx2 but not otx2
% staining)
include = inputdlg('Include OTX2?'); include = str2double(include);
if include == 1
    otx2data(:,i) = alldata(:,10); otx2stdev(:,i) = alldata(:,11);
elseif include == 0
    otx2data(:,i) = NaN; otx2stdev(:,i) = NaN;
end

end

% average cell densities (take avg of each row)
avgdensity = mean(density,2,'omitnan');
stdevdensity = std(density,0,2,'omitnan'); % 'w' = 0 normalizes by n-1 (default)
semdensity = stdevdensity/sqrt(n);

% take weighted average of intensities
totalcounts = sum(count,2);
weightedcount = count./totalcounts;

% find average (weighted) intensities
avghoeschtdata = sum(hoeschtdata.*weightedcount,2,'omitnan');
avgfoxa2data = sum(foxa2data.*weightedcount,2,'omitnan');
avgcdx2data = sum(cdx2data.*weightedcount,2,'omitnan');
avgotx2data = sum(otx2data.*weightedcount,2,'omitnan');

% propagate uncertainty using pooled variance
pooledweight = (count - 1)./(totalcounts - n);

avghoeschtstdev = (sum((hoeschtstdev.^2).*pooledweight,2,'omitnan')).^0.5;
avgfoxa2stdev = (sum((foxa2stdev.^2).*pooledweight,2,'omitnan')).^0.5;
avgcdx2stdev = (sum((cdx2stdev.^2).*pooledweight,2,'omitnan')).^0.5;
avgotx2stdev = (sum((otx2stdev.^2).*pooledweight,2,'omitnan')).^0.5;

% calculate sem for plotting
avghoeschtsem = avghoeschtstdev/sqrt(n);
avgfoxa2sem = avgfoxa2stdev/sqrt(n);
avgcdx2sem = avgcdx2stdev/sqrt(n);
avgotx2sem = avgotx2stdev/sqrt(n);

%% Generate Plots

% choose save destination
folder_name = uigetdir(); cd(folder_name);

figure; % plot density
shadedErrorBar(x, avgdensity, semdensity,'lineProps','-b','transparent',1);
title('Nuclear density'); xlabel('x (\mum)'); ylabel('Nuclei per mm^2')
xlim([0 1000]); ylim([0 2500]) 
saveas(gcf,'Average Nuclear Density','epsc') % save in vector format
saveas(gcf,'Average Nuclear Density','tiffn') % save in tif format

figure; % plot Hoescht
shadedErrorBar(x, avghoeschtdata/max(avghoeschtdata),...
    avghoeschtsem/max(avghoeschtdata),'lineProps','-b','transparent',1);
title('Normalized Hoescht Intensity')
xlabel('x (\mum)'); ylabel('Normalized Hoescht Intensity'); 
xlim([0 1000]); ylim([0 1.5])
% save in vector and tif format
saveas(gcf,'Avg Norm Hoescht Intensity','epsc'); 
saveas(gcf,'Avg Norm Hoescht Intensity','tiffn');

figure; % plot foxa2
shadedErrorBar(x, avgfoxa2data/max(avgfoxa2data),...
    avgfoxa2sem/max(avgfoxa2data),'lineProps','-g','transparent',1);
title('Normalized FOXA2 Intensity')
xlabel('x (\mum)'); ylabel('Normalized FOXA2 Intensity'); 
xlim([0 1000]); ylim([0 1.5])
% save in vector and tif format
saveas(gcf,'Avg Norm FOXA2 Intensity','epsc'); 
saveas(gcf,'Avg Norm FOXA2 Intensity','tiffn');

figure; % plot cdx2
shadedErrorBar(x, avgcdx2data/max(avgcdx2data),...
    avgcdx2sem/max(avgcdx2data),'lineProps','-m','transparent',1);
title('Normalized CDX2 Intensity')
xlabel('x (\mum)'); ylabel('Normalized CDX2 Intensity'); 
xlim([0 1000]); ylim([0 1.5])
% save in vector and tif format
saveas(gcf,'Avg Norm CDX2 Intensity','epsc'); 
saveas(gcf,'Avg Norm CDX2 Intensity','tiffn');

figure; % plot otx2
shadedErrorBar(x, avgotx2data/max(avgotx2data),...
    avgotx2sem/max(avgotx2data),'lineProps','-y','transparent',1);
title('Normalized OTX2 Intensity')
xlabel('x (\mum)'); ylabel('Normalized OTX2 Intensity'); 
xlim([0 1000]); ylim([0 1.5])
% save in vector and tif format
saveas(gcf,'Avg Norm OTX2 Intensity','epsc'); 
saveas(gcf,'Avg Norm OTX2 Intensity','tiffn');

%% Save Data

output.x = x;
output.avgdensity = avgdensity;
output.semdensity = semdensity;
output.hoeschtint = avghoeschtdata;
output.hoeschtsem = avghoeschtsem;
output.foxa2int = avgfoxa2data;
output.foxa2sem = avgfoxa2sem;
output.cdx2int = avgcdx2data;
output.cdx2sem = avgcdx2sem;
output.otx2int = avgotx2data;
output.otx2sem = avgotx2sem;

fileID = fopen('AvgAnalysis.txt','w');
writetable(struct2table(output),'AvgAnalysis.txt')
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
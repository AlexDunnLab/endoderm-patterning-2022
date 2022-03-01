function [likelihood_func_CDX2, likelihood_func_OTX2] = onegene_interp_gut_data_gaussians(points, nbins)
%% Return the Gaussian approximation for the distribution of expression
% levels over positions, for each gene separately.
% INPUT: 
% * points, Table with columns: x, CDX2, OTX2
% * nbins, int, number of bins along the x axis.
% OUTPUT:
% * likelihood_funcs, function handle for a function that takes 2 arrays 
%                     x, cdx2 (or otx2) and returns the likelihood.

%% Bin points into the designated number of bins
[bins,binedges] = discretize(points.x,nbins);

%% means and Covariance matrices at every bin position
Cs_OTX2 = zeros(1,nbins);
Cs_CDX2 = zeros(1,nbins);
mus_OTX2  = zeros(1,nbins);
mus_CDX2  = zeros(1,nbins);
bincenters = (binedges(1:end-1)+binedges(2:end))./2;

for i=1:nbins
    thisbin_points = points(bins==i,2:3);
    Cs_OTX2(i) = var([thisbin_points.OTX2]); %covariance 
    Cs_CDX2(i) = var([thisbin_points.CDX2]); %covariance 
    mus_OTX2(i) = mean(thisbin_points.OTX2); %mean
    mus_CDX2(i) = mean(thisbin_points.CDX2); %mean
end

%% return a function that interpolates the means and covariances 
% along the interval

    function y2 = gut_gaussian_approx_CDX2(x,cdx2,otx2)
        %gaussian approximation to the joint distribution of (x,cdx2)

        y2 = zeros(length(x),1);
        for idx = 1:length(x)
            C = interp1([-1,bincenters,2],[Cs_CDX2(1), Cs_CDX2 ,Cs_CDX2(end)],x(idx));
            mu = interp1([-1,bincenters, 2],[mus_CDX2(1), mus_CDX2 ,mus_CDX2(end)],x(idx));
            y2(idx) = normpdf(cdx2(idx),squeeze(mu),sqrt(squeeze(C)));
        end
    end

    function y2 = gut_gaussian_approx_OTX2(x,otx2,cdx2)
            %gaussian approximation to the joint distribution of (x,otx2)

            y2 = zeros(length(x),1);
            for idx = 1:length(x)
                C = interp1([-1,bincenters,2],[Cs_OTX2(1), Cs_OTX2 ,Cs_OTX2(end)],x(idx));
                mu = interp1([-1,bincenters, 2],[mus_OTX2(1), mus_OTX2 ,mus_OTX2(end)],x(idx));
                y2(idx) = normpdf(otx2(idx),squeeze(mu),sqrt(squeeze(C)));
            end
        end

    likelihood_func_CDX2 = @gut_gaussian_approx_CDX2;
    likelihood_func_OTX2 = @gut_gaussian_approx_OTX2;
    end
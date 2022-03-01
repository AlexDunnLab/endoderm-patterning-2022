function likelihood_func = interp_gut_data_gaussians(points, nbins)
%% Return the Gaussian approximation for the distribution of expression
% levels over positions.
% INPUT: 
% * points, Table with columns: x, CDX2, OTX2
% * nbins, int, number of bins along the x axis.
% OUTPUT:
% * likelihood_func, function handle for a function that takes 3 arrays 
%                     x, cdx2, otx2 and returns the likelihood.

%% Bin points into the designated number of bins
[bins,binedges] = discretize(points.x,nbins);

%% means and Covariance matrices at every bin position
Cs = zeros(nbins,2,2);
mus  = zeros(nbins,2);
bincenters = (binedges(1:end-1)+binedges(2:end))./2;

for i=1:nbins
    thisbin_points = points(bins==i,2:3);
    Cs(i,:,:) = cov([thisbin_points.CDX2 thisbin_points.OTX2]); %covariance matrix 
    mus(i,:) = mean([thisbin_points.CDX2 thisbin_points.OTX2]); %mean
end

%% return a function that interpolates the means and covariances 
% along the interval

    function y2 = gut_gaussian_approx(x,cdx2,otx2)
        %gaussian approximation to the joint distribution of (x,cdx2,otx2)

        y2 = zeros(length(x),1);
        for idx = 1:length(x)
            C = interp1([-1,bincenters,2],cat(1,Cs(1,:,:), Cs ,Cs(end,:,:)),x(idx));
            mu = interp1([-1,bincenters, 2],cat(1,mus(1,:), mus ,mus(end,:)),x(idx));
            y2(idx) = mvnpdf([cdx2(idx),otx2(idx)],squeeze(mu),squeeze(C));
        end
    end

likelihood_func = @gut_gaussian_approx;
end
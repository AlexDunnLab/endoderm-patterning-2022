function likelihood_func = interp_gut_data_bernoulli(points,nbins,CDX2thresh,OTX2thresh)
%% Return the Bernoulli approximation for the distribution of expression
% levels over positions, jointly for the genes.
% INPUT: 
% * points, Table with columns: x, CDX2, OTX2
% * nbins, int, number of bins along the x axis.
% OUTPUT:
% * likelihood_funcs, function handle for a function that takes 2 arrays 
%                     x, cdx2 (or otx2) and returns the likelihood.

%     CDX2thresh = 0.1;graythresh(points.CDX2);
%     OTX2thresh = 0.35;graythresh(points.OTX2);
       
    [bins,binedges] = discretize(points.x,nbins);
    bincenters = (binedges(1:end-1)+binedges(2:end))./2;
    probmat = zeros(length(bincenters),2,2);
    for i=1:nbins
        thisbin_points = points(bins==i,1:3);
        npoints = length(thisbin_points.x);
        otx2high = thisbin_points.OTX2>OTX2thresh;
        cdx2high = thisbin_points.CDX2>CDX2thresh;
        probmat(i,1,1) = sum(otx2high.*cdx2high)/npoints;
        probmat(i,1,2) = sum((1-otx2high).*cdx2high)/npoints;
        probmat(i,2,1) = sum((otx2high).*(1-cdx2high))/npoints;
        probmat(i,2,2) = 1-probmat(i,1,1)-probmat(i,1,2)-probmat(i,2,1);
    end
    
    %% return a function that interpolates the means and covariances 
    % along the interval

    function y2 = gut_bernoulli_approx(x,cdx2,otx2)
        %bernoulli approximation to the joint distribution of (x,cdx2)

        y2 = zeros(length(x),1);
        for idx = 1:length(x)
            p = squeeze(interp1([-1,bincenters,2],cat(1,probmat(1,:,:), probmat ,probmat(end,:,:)),x(idx)));
            sssss = sum(p);
            if (cdx2(idx)>CDX2thresh)
                if (otx2(idx)>OTX2thresh)
                    y2(idx) = p(1,1);
                else
                    y2(idx) = p(1,2);
                end
            else
                if (otx2(idx)>OTX2thresh)
                    y2(idx) = p(2,1);
                else
                    y2(idx) = p(2,2);
                end
            end
            
        end
    end


    likelihood_func = @gut_bernoulli_approx;

end
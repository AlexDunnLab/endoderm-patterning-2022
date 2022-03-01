function [likelihood_func_CDX2, likelihood_func_OTX2] = onegene_interp_gut_data_bernoulli(points,nbins)
%% Return the Bernoulli approximation for the distribution of expression
% levels over positions, for each gene separately.
% INPUT: 
% * points, Table with columns: x, CDX2, OTX2
% * nbins, int, number of bins along the x axis.
% OUTPUT:
% * likelihood_funcs, function handle for a function that takes 2 arrays 
%                     x, cdx2 (or otx2) and returns the likelihood.

    CDX2thresh = 0.1518;%max(points.CDX2).*(graythresh(points.CDX2/max(points.CDX2)))
    OTX2thresh = 0.3389;%max(points.OTX2).*(graythresh(points.OTX2/max(points.OTX2)))
    
    [bins,binedges] = discretize(points.x,nbins);
    bincenters = (binedges(1:end-1)+binedges(2:end))./2;
    
    for i=1:nbins
        thisbin_points = points(bins==i,1:3);
        Ps_OTX2(i) = sum([thisbin_points.OTX2>OTX2thresh])./length(thisbin_points.x); %fraction
        Ps_CDX2(i) = sum([thisbin_points.CDX2>CDX2thresh])./length(thisbin_points.x); %fraction 
    end
    
    %% return a function that interpolates the means and covariances 
    % along the interval

    function y2 = gut_bernoulli_approx_CDX2(x,cdx2)
        %bernoulli approximation to the joint distribution of (x,cdx2)

        y2 = zeros(length(x),1);
        for idx = 1:length(x)
            p = interp1([-1,bincenters,2],[Ps_CDX2(1), Ps_CDX2 ,Ps_CDX2(end)],x(idx));
            if (cdx2>CDX2thresh)
                y2(idx) = p;
            else
                y2(idx) = 1-p;
            end
            
        end
    end

    function y2 = gut_bernoulli_approx_OTX2(x,otx2)
            %gaussian approximation to the joint distribution of (x,otx2)

            y2 = zeros(length(x),1);
            for idx = 1:length(x)
                p = interp1([-1,bincenters,2],[Ps_OTX2(1), Ps_OTX2 ,Ps_OTX2(end)],x(idx));
                if (otx2>OTX2thresh)
                    y2(idx) = p;
                else
                    y2(idx) = 1-p;
                end
            
            end
    end

    likelihood_func_CDX2 = @gut_bernoulli_approx_CDX2;
    likelihood_func_OTX2 = @gut_bernoulli_approx_OTX2;

end
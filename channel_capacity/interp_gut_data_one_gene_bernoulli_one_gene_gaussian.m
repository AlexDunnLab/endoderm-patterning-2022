function likelihood_func = interp_gut_data_one_gene_bernoulli_one_gene_gaussian(points,nbins,CDX2thresh)
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
    probmat = zeros(length(bincenters),1); %probabilities that CDX2>thresh
    
    Cs_high = zeros(nbins,1);
    mus_high  = zeros(nbins,1);
    Cs_low = zeros(nbins,1);
    mus_low  = zeros(nbins,1);
    for i=1:nbins
        thisbin_points = points(bins==i,1:3);
        npoints = length(thisbin_points.x);
        cdx2high = thisbin_points.CDX2>CDX2thresh;
        Cs_high(i) = var(thisbin_points.OTX2(thisbin_points.CDX2>CDX2thresh)); %variances 
        mus_high(i) = mean(thisbin_points.OTX2(thisbin_points.CDX2>CDX2thresh)); %means
        Cs_low(i) = var(thisbin_points.OTX2(thisbin_points.CDX2<CDX2thresh)); %variances 
        mus_low(i) = mean(thisbin_points.OTX2(thisbin_points.CDX2<CDX2thresh)); %means
        probmat(i) = sum(cdx2high)/npoints;
    end


    %% return a function that interpolates the means and covariances 
    % along the interval
    
    function y2 = gut_bernoulligaussian_approx(x,cdx2,otx2)
        %bernoulli approximation to the joint distribution of (x,cdx2)
        
        y2 = zeros(length(x),1);
        for idx = 1:length(x)
            p = squeeze(interp1([-1,bincenters,2],cat(1,probmat(1), probmat ,probmat(end)),x(idx)));
            sssss = sum(p);
            if (cdx2(idx)>CDX2thresh)
                prob = p;
                C = interp1([-1,bincenters,2],cat(1,Cs_high(1), Cs_high ,Cs_high(end)),x(idx));
                mu = interp1([-1,bincenters, 2],cat(1,mus_high(1), mus_high ,mus_high(end)),x(idx));
            else
                prob = 1-p;
                C = interp1([-1,bincenters,2],cat(1,Cs_low(1), Cs_low ,Cs_low(end)),x(idx));
                mu = interp1([-1,bincenters, 2],cat(1,mus_low(1), mus_low ,mus_low(end)),x(idx));
            end
            y2(idx) = prob .* normpdf(otx2(idx),squeeze(mu),sqrt(squeeze(C)));
        end
    end


    likelihood_func = @gut_bernoulligaussian_approx;

end
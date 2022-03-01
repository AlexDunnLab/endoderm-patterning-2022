function decoder_map = optimalDecoder(points,nbins)
%% Given the profile
% return a heatmap of p(x*|x) = 
% = int(dotx2 dcdx2 p(x*|(otx2,cdx2))*p((OTX2,CDX2)|x)) 
% = int(dotx2 dcdx2 p((otx2,cdx2)|x*)*p(x*)/p(otx2,cdx2)*p((OTX2,CDX2)|x))
% \propto int(dotx2 dcdx2 p((otx2,cdx2)|x*) p((otx2,cdx2)|x) ) 
% as a function of x

%% Generate the whole joint distribution
% func = interp_gut_data_gaussians(points, nbins);
% func = interp_gut_data_bernoulli(points, nbins,0.2,0.4);
% func = interp_gut_data_one_gene_bernoulli_one_gene_gaussian(points, nbins,0.2);
[~,func] = onegene_interp_gut_data_gaussians(points,nbins);
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

%% Generate decoder map
decoder_map = zeros(nbins);

dc = 0.05;

for i = 1:length(bincenters)
    [cdx2grid,otx2grid] = meshgrid(-0.1:dc:1.1,-0.1:dc:1.1);
    xstargrid = bincenters(i)+0.*cdx2grid;
    p_xstaroc = func(xstargrid(:),cdx2grid(:),otx2grid(:));
    p_oc = 0*func(xstargrid(:),cdx2grid(:),otx2grid(:));

%     heatmap(reshape(p_oc,size(cdx2grid)));
%     figure;
    for j = 1:length(bincenters)
        xgrid = bincenters(j)+0.*cdx2grid;
        p_oc = p_oc + 1/nbins * func(xgrid(:),cdx2grid(:),otx2grid(:));
    end
    for j = 1:length(bincenters)
        xgrid = bincenters(j)+0.*cdx2grid;
%         xgrid = 0.5*bincenters(j)+0.*cdx2grid;
        p_xoc = func(xgrid(:),cdx2grid(:),otx2grid(:));
        decoder_map(i,j) = nansum(nansum(p_xoc.*p_xstaroc./p_oc));
    end
end

end
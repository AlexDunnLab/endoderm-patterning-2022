%% Fake Data

x = linspace(0,1,2000)';
sharpness = 5;
otx2 = heaviside(linspace(sharpness,-sharpness,2000))'/3+0.5+0.2.*randn(2000,1);
cdx2 = heaviside(linspace(-sharpness,sharpness,2000))'/3+0.5+0.2.*randn(2000,1);
data = table(x,otx2,cdx2);
data.Properties.VariableNames = {'x','OTX2','CDX2'};


% %% Generate gaussian fits
% func = interp_gut_data_gaussians(data,20);
% [func_CDX2, func_OTX2] = onegene_interp_gut_data_gaussians(data,20);


%% see the raw data for OTX2, CDX2
% scatter(data.x,data.OTX2)
% hold on;
% scatter(data.x,data.CDX2)

% 
%% Compute channel capacity

[funcCDX, funcOTX ] = onegene_interp_gut_data_gaussians(data,20);

dx = 0.02;
dcdx = 0.1;
xpoints = 0:dx:1;
cdxpoints = -0.1:dcdx:1.1;
full_joint_dist = zeros(length(xpoints),length(cdxpoints).^2);
for i=1:length(xpoints)
    [xmesh,CDX2mesh, OTX2mesh] = meshgrid(xpoints(i),cdxpoints,cdxpoints);
    funcmesh = funcCDX(xmesh(:),CDX2mesh(:),OTX2mesh(:))';
    full_joint_dist(i,:) = funcmesh./sum(funcmesh);
end
imagesc(full_joint_dist)

[Cap,r] = BlahutArimoto(full_joint_dist);
disp('CDX2 only:');disp(Cap);

func = interp_gut_data_bernoulli(data,20,CDXthresh,OTXthresh);

dx = 0.02;
dcdx = 0.1;
xpoints = 0:dx:1;
cdxpoints = -0.1:dcdx:1.1;
full_joint_dist = zeros(length(xpoints),length(cdxpoints).^2);
for i=1:length(xpoints)
    [xmesh,CDX2mesh, OTX2mesh] = meshgrid(xpoints(i),cdxpoints,cdxpoints);
    funcmesh = funcOTX(xmesh(:),CDX2mesh(:),OTX2mesh(:))';
    full_joint_dist(i,:) = funcmesh./sum(funcmesh);
end
imagesc(full_joint_dist)

[Cap,r] = BlahutArimoto(full_joint_dist);
disp('OTX2 only:');disp(Cap);


% CDXthresh = xf(1);
% OTXthresh = xf(2);
% CDXthresh = .2;
% OTXthresh = .5;
% func = interp_gut_data_bernoulli(data,20,CDXthresh,OTXthresh);
% 
% dx = 0.02;
% dcdx = 0.1;
% xpoints = 0:dx:1;
% cdxpoints = -0.1:dcdx:1.1;
% full_joint_dist = zeros(length(xpoints),length(cdxpoints).^2);
% for i=1:length(xpoints)
%     [xmesh,CDX2mesh, OTX2mesh] = meshgrid(xpoints(i),cdxpoints,cdxpoints);
%     funcmesh = func(xmesh(:),CDX2mesh(:),OTX2mesh(:))';
%     full_joint_dist(i,:) = funcmesh./sum(funcmesh);
% end
% imagesc(full_joint_dist)
% 
% [Cap,r] = BlahutArimoto(full_joint_dist);
% disp('both binary:');disp(Cap);
% 
% full_joint_dist = zeros(length(xpoints),length(cdxpoints).^2);
% func = interp_gut_data_one_gene_bernoulli_one_gene_gaussian(data,20,CDXthresh);
% for i=1:length(xpoints)
%     [xmesh,CDX2mesh, OTX2mesh] = meshgrid(xpoints(i),cdxpoints,cdxpoints);
%     funcmesh = func(xmesh(:),CDX2mesh(:),OTX2mesh(:))';
%     full_joint_dist(i,:) = funcmesh./nansum(funcmesh);
% end
% imagesc(full_joint_dist)
% full_joint_dist(isnan(full_joint_dist))=0;
% [Cap,r] = BlahutArimoto(full_joint_dist);
% disp('CDX2 binary:');disp(Cap);
% 
% data2 = data(:,:);
% data2.OTX2 = data.CDX2;
% data2.CDX2 = data.OTX2;
% 
% func = interp_gut_data_one_gene_bernoulli_one_gene_gaussian(data2,20,OTXthresh);
% full_joint_dist = zeros(length(xpoints),length(cdxpoints).^2);
% for i=1:length(xpoints)
%     [xmesh,CDX2mesh, OTX2mesh] = meshgrid(xpoints(i),cdxpoints,cdxpoints);
%     funcmesh = func(xmesh(:),CDX2mesh(:),OTX2mesh(:))';
%     full_joint_dist(i,:) = funcmesh./sum(funcmesh);
% end
% imagesc(full_joint_dist)
% 
% [Cap,r] = BlahutArimoto(full_joint_dist);
% disp('OTX2 binary:');disp(Cap);

%% Test that the marginal distribution of CDX2 makes sense.
% figure;
% dx = 0.01;
% dcdx = 0.01;
% xpoints = 0:dx:1;
% cdxpoints = -0.1:dcdx:1.1;
% [xtest,cdx2test] = meshgrid([0],cdxpoints);
% hmapdata = func_OTX2(xtest(:),cdx2test(:));
% for i = xpoints(2:end)
%     [xtest,cdx2test] = meshgrid([i],cdxpoints);
%     hmapdata = hmapdata + func_OTX2(xtest(:),cdx2test(:));
% end
% plot(cdxpoints,squeeze((hmapdata)));
% 
% figure;
% [xtest,cdx2test] = meshgrid(xpoints,cdxpoints);
% surf(cdx2test,xtest,reshape(func_OTX2(xtest(:),cdx2test(:)), length(cdxpoints),length(xpoints)));
% totalint = sum(sum(sum(hmapdata)))*dx*dcdx


% Compute a decoder map

% data2 = data(:,:);
% data2.OTX2 = data.CDX2;
% data2.CDX2 = data.OTX2;
% decoder_map = optimalDecoder(data2,50);
% 
% %normalize decoder map:
% decoder_map = decoder_map./(sum(decoder_map));

% figure;
% colormap('gray');
% imagesc(imcomplement(decoder_map));

% % 
% Compute a mutual information
MI_resolution = 0.03;
MI_CDX2alone = onegene_mutual_information_riemann(func_CDX2,MI_resolution,-1,2,MI_resolution)
MI_OTX2alone = onegene_mutual_information_riemann(func_OTX2,MI_resolution,-1,2,MI_resolution)
MI_both = mutual_information_riemann(func,MI_resolution,0,1,MI_resolution,-0,1,MI_resolution)
MI_withprior = mutual_information_riemann_hill_gradient(func,MI_resolution,0,1,MI_resolution,0,1,MI_resolution,testpdf)
figure;
% 
% x = linspace(0.1,0.9,100);
% [xgrid,cdx2grid,otx2grid ]= meshgrid(x,x,[0.7]);
% 
% hmapdata = (reshape(func(xgrid(:),cdx2grid(:),otx2grid(:)),100,100));
% hmapdata2 = (reshape(func(xgrid(:),otx2grid(:),cdx2grid(:)),100,100));
% heatmap(hmapdata./sum(hmapdata));
% figure;
% heatmap(hmapdata2./sum(hmapdata2));


% %% Generate Bernoulli fits
% 
% [func_CDX2, func_OTX2] = onegene_interp_gut_data_bernoulli(data,200);
% MI_resolution = 0.001;
% MI_CDX2alone = onegene_mutual_information_bernoulli(func_CDX2,MI_resolution)
% MI_OTX2alone = onegene_mutual_information_bernoulli(func_OTX2,MI_resolution)
% 
% opt_thresh = returnObj(data,MI_resolution);
% % xf = fminsearch(opt_thresh,[0.1,0.35])
% CDXthresh = xf(1);
% OTXthresh = xf(2);
% func = interp_gut_data_bernoulli(data,20,CDXthresh,OTXthresh);


% MI_both = mutual_information_bernoulli(func,MI_resolution)
% %% Test that the marginal distribution of OTX2 x CDX2 makes sense.
% figure;
% dx = 0.1;
% dcdx = 0.05;
% dotx=0.05;
% xpoints = 0:dx:1;
% cdxpoints = -0.3:dcdx:1.3;
% otxpoints = -0.3:dotx:1.3;
% [xtest,cdx2test,otx2test] = meshgrid([0.1],cdxpoints,otxpoints);
% hmapdata = reshape(func(xtest(:),cdx2test(:),otx2test(:)),length(cdxpoints),length(otxpoints));
% for i = xpoints(2:end)
%     [xtest,cdx2test,otx2test] = meshgrid([i],cdxpoints,otxpoints);
%     hmapdata = hmapdata + reshape(func(xtest(:),cdx2test(:),otx2test(:)),length(cdxpoints),length(otxpoints));
% end
% surf(cdxpoints,otxpoints,squeeze((hmapdata)));
% totalint = sum(sum(sum(hmapdata)))*dx*dotx*dcdx
% function Y = returnObj(data1,MI_resolution)
%     function MI_both = optimize_threshes(x1)
%         CDXthresh = x1(1);
%         OTXthresh = x1(2);
%         func = interp_gut_data_bernoulli(data1,20,CDXthresh,OTXthresh);
%         MI_both = -mutual_information_bernoulli(func,MI_resolution)
%     end
%     Y = @optimize_threshes;
% end
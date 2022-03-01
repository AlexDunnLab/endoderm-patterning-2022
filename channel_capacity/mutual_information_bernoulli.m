function MI = mutual_information_bernoulli(joint_pdf,dx)
    x_points = (0:dx:1)';
    n_x = length(x_points);
    
    
    %% marginal probability mass function of two genes being high/low
    [xgrid_marg,cdx2_marg, otx2_marg] = meshgrid(x_points(1),[-1,1],[-1,1]);
    mpmf = dx*joint_pdf(xgrid_marg(:),cdx2_marg(:),otx2_marg(:));
    for idx = x_points(2:end)'
        [xgrid_marg,cdx2_marg, otx2_marg] = meshgrid(idx,[-1,1],[-1,1]);
        new = joint_pdf(xgrid_marg(:),cdx2_marg(:),otx2_marg(:));
        mpmf  = mpmf + dx*new;
    end
    
    %% Compute mutual information
    [xgrid_marg,cdx2_marg, otx2_marg] = meshgrid(x_points(1),[0,1],[0,1]);
    jpdf = joint_pdf(xgrid_marg(:),cdx2_marg(:),otx2_marg(:));
    MI = dx.*jpdf.*log2(jpdf./mpmf);
    MI(isnan(MI))=0;
    for idx = x_points(2:end)'
        [xgrid_marg,cdx2_marg, otx2_marg] = meshgrid(idx,[0,1],[0,1]);
        jpdf = joint_pdf(xgrid_marg(:),cdx2_marg(:),otx2_marg(:));
        new = dx*jpdf.*log2(jpdf./mpmf);
        new(isnan(new))=0;
        MI = MI+new;
    end    
    
%     
%     mpdf = mpdf./sum(sum(sum(mpdf)));
%     mpdfres = reshape(mpdf,n_x,n_cdx2,n_otx2);
%     heatmap(cdx2_points,otx2_points,squeeze(mpdfres(1,:,:)));
    
%     jpdf_sanity = dx*dcdx2*dotx2*sum(sum(sum(jpdf)))
%     mpdf_sanity = dg*sum(sum(sum(mprob)))
    
%     matforint = dx*dotx2*dcdx2*jpdf.*(log2(jpdf./repmat(mpdf,length(x_points),1)));
    
    MI = nansum(MI);
    
end
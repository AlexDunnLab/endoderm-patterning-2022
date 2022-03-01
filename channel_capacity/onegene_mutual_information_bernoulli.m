function MI = onegene_mutual_information_bernoulli(joint_pdf,dx)
    x_points = (0:dx:1)';
    n_x = length(x_points);
    
    
    %% total probability of gene being high
    [xgrid_marg,genegrid_marg] = meshgrid(x_points(:),1);
    mprob = dx*sum(sum(joint_pdf(xgrid_marg(:),genegrid_marg(:))));
    
    %% Compute mutual information
    [xgrid_marg,genegrid_marg] = meshgrid(x_points(:),1);
    jpdf = joint_pdf(xgrid_marg(:),genegrid_marg(:));
    MI = dx.*nansum(nansum((jpdf.*log2(jpdf./mprob) + (ones('like', jpdf)-jpdf).*log2((ones('like', jpdf)-jpdf)./(1-mprob)))));
%     jpdf = joint_pdf(xgrid(:),cdx2grid(:),otx2grid(:));
%     jpdf = jpdf./sum(sum(sum(jpdf)));

    
    
%     
%     mpdf = mpdf./sum(sum(sum(mpdf)));
%     mpdfres = reshape(mpdf,n_x,n_cdx2,n_otx2);
%     heatmap(cdx2_points,otx2_points,squeeze(mpdfres(1,:,:)));
    
%     jpdf_sanity = dx*dcdx2*dotx2*sum(sum(sum(jpdf)))
%     mpdf_sanity = dg*sum(sum(sum(mprob)))
    
%     matforint = dx*dotx2*dcdx2*jpdf.*(log2(jpdf./repmat(mpdf,length(x_points),1)));
    
    MI = nansum(MI);
    
end
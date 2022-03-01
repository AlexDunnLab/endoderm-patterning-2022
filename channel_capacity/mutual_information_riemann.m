function MI = mutual_information_riemann(joint_pdf,dx,cdx2_min,cdx2_max,dcdx2,otx2_min,otx2_max,dotx2)
    x_points = (0:dx:1)';
    n_x = length(x_points)
    
    cdx2_points = (cdx2_min:dcdx2:cdx2_max)';
    n_cdx2 = length(cdx2_points);
    
    otx2_points = (otx2_min:dotx2:otx2_max)';
    n_otx2 = length(otx2_points);
    
    %% Compute marginal distribution of expression
    [xgrid_marg,cdx2grid_marg,otx2grid_marg] = meshgrid([x_points(1)],cdx2_points,otx2_points);
    mpdf = dx*joint_pdf(xgrid_marg(:),cdx2grid_marg(:),otx2grid_marg(:));
    for idx = x_points(2:end)'
        [xgrid_marg,cdx2grid_marg,otx2grid_marg] = meshgrid([idx],cdx2_points,otx2_points);
        new = joint_pdf(xgrid_marg(:),cdx2grid_marg(:),otx2grid_marg(:));
        mpdf = mpdf + dx*new;
    end
    
    %% Compute mutual information
    [xgrid_marg,cdx2grid_marg,otx2grid_marg] = meshgrid([x_points(1)],cdx2_points,otx2_points);
    jpdf = joint_pdf(xgrid_marg(:),cdx2grid_marg(:),otx2grid_marg(:));
    MI = dcdx2*dotx2*dx.*jpdf.*log2(jpdf./mpdf);
    for idx = x_points(2:end)'
        [xgrid_marg,cdx2grid_marg,otx2grid_marg] = meshgrid([idx],cdx2_points,otx2_points);
        jpdf = joint_pdf(xgrid_marg(:),cdx2grid_marg(:),otx2grid_marg(:));
        MI = MI + dcdx2*dotx2*dx*  jpdf   .*  log2(jpdf./mpdf);
    end
%     jpdf = joint_pdf(xgrid(:),cdx2grid(:),otx2grid(:));
%     jpdf = jpdf./sum(sum(sum(jpdf)));

    
    
%     
%     mpdf = mpdf./sum(sum(sum(mpdf)));
%     mpdfres = reshape(mpdf,n_x,n_cdx2,n_otx2);
%     heatmap(cdx2_points,otx2_points,squeeze(mpdfres(1,:,:)));
    
%     jpdf_sanity = dx*dcdx2*dotx2*sum(sum(sum(jpdf)))
    mpdf_sanity = dcdx2*dotx2*sum(sum(sum(mpdf)))
    
%     matforint = dx*dotx2*dcdx2*jpdf.*(log2(jpdf./repmat(mpdf,length(x_points),1)));
    
    MI = nansum(MI);
    
end
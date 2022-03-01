function MI = onegene_mutual_information_riemann(joint_pdf,dx,gene_min,gene_max,dg)
    x_points = (0:dx:1)';
    
    %step function of xes
%     x_points = cat(2,zeros(1,50),ones(1,50))';
    n_x = length(x_points)
    
    gene_points = (gene_min:dg:gene_max)';
    n_gene = length(gene_points);
    
    
    %% Compute marginal distribution of expression for each gene separately
    [xgrid_marg,genegrid_marg] = meshgrid([x_points(1)],gene_points(:));
    mpdf = dx*joint_pdf(xgrid_marg(:),genegrid_marg(:));
    for idx = x_points(2:end)'
        [xgrid_marg,genegrid_marg] = meshgrid([idx],gene_points(:));
        new = joint_pdf(xgrid_marg(:),genegrid_marg(:));
        mpdf = mpdf + dx*new;
    end
    
    %% Compute mutual information
    [xgrid_marg,genegrid_marg] = meshgrid([x_points(1)],gene_points(:));
    jpdf = joint_pdf(xgrid_marg(:),genegrid_marg(:));
    MI = dg*dx.*jpdf.*log2(jpdf./mpdf);
    for idx = x_points(2:end)'
        [xgrid_marg,genegrid_marg] = meshgrid([idx],gene_points(:));
        jpdf = joint_pdf(xgrid_marg(:),genegrid_marg(:));
        MI = MI + dg*dx*  jpdf   .*  log2(jpdf./mpdf);
    end
%     jpdf = joint_pdf(xgrid(:),cdx2grid(:),otx2grid(:));
%     jpdf = jpdf./sum(sum(sum(jpdf)));

    
    
%     
%     mpdf = mpdf./sum(sum(sum(mpdf)));
%     mpdfres = reshape(mpdf,n_x,n_cdx2,n_otx2);
%     heatmap(cdx2_points,otx2_points,squeeze(mpdfres(1,:,:)));
    
%     jpdf_sanity = dx*dcdx2*dotx2*sum(sum(sum(jpdf)))
    mpdf_sanity = dg*sum(sum(sum(mpdf)))
    
%     matforint = dx*dotx2*dcdx2*jpdf.*(log2(jpdf./repmat(mpdf,length(x_points),1)));
    
    MI = nansum(MI);
    
end
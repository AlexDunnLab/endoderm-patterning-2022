function MI = mutual_information_riemann_hill_gradient(joint_pdf,dx,cdx2_min,cdx2_max,dcdx2,otx2_min,otx2_max,dotx2,x_prior_dist)
    x_points = (0:dx:1)';
    n_x = length(x_points)
    
    cdx2_points = (cdx2_min:dcdx2:cdx2_max)';
    n_cdx2 = length(cdx2_points);
    
    otx2_points = (otx2_min:dotx2:otx2_max)';
    n_otx2 = length(otx2_points);
    
    xprior_norm = sum(x_prior_dist(x_points));
    
    conditional_prob_mat = zeros(n_x,n_cdx2,n_otx2);
    
    xprior = x_prior_dist(x_points)/xprior_norm;
    
    for i=1:n_x
        for c = 1:n_cdx2
            for o = 1:n_otx2
                conditional_prob_mat(i,c,o) = joint_pdf(x_points(i),cdx2_points(c),otx2_points(o));
            end
        end
    end
    
    % normalize the x bins
    
    for i = 1:n_x
        conditional_prob_mat(i,:,:) = conditional_prob_mat(i,:,:)/sum(sum(conditional_prob_mat(i,:,:)));
    end
    
    mpdf = 0*conditional_prob_mat(1,:,:);
    
    for i = 1:n_x
        mpdf = mpdf + xprior(i)*conditional_prob_mat(i,:,:);
    end
    
    MI_to_sum = 0*conditional_prob_mat;
    for i = 1:n_x
       MI_to_sum(i,:,:) = conditional_prob_mat(i,:,:) * xprior(i).*log2(conditional_prob_mat(i,:,:)./mpdf);
    end
    
    MI = nansum(nansum(nansum(MI_to_sum)));
    
end
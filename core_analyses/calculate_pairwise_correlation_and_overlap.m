function imdata = calculate_pairwise_correlation_and_overlap(Iall,Isyn,Isyn_trans,ctrlshift_xy)
%Calculates 3D and 2D correlation between synaptic channels
%Control image is the image shifted to check for spurious correlation

correlations = zeros(size(Iall,2),size(Iall,2));
overlaps = zeros(size(Iall,2),size(Iall,2));
correlations_ctrl = zeros(size(Iall,2),size(Iall,2));
overlaps_ctrl = zeros(size(Iall,2),size(Iall,2));
correlations_2D = zeros(size(Iall,2),size(Iall,2));
correlations_2D_ctrl = zeros(size(Iall,2),size(Iall,2));

for c1 = 1:size(Iall,2)
    for c2 = c1+1:size(Iall,2)
        %3D correlation
        correlations(c1,c2) = corr(Isyn{c1}(:),Isyn{c2}(:));
        intersection = min(Isyn{c1},Isyn{c2});
        union = max(Isyn{c1},Isyn{c2});
        overlaps(c1,c2) = sum(intersection(:))./sum(union(:));
        
        translated = Isyn_trans{c2}(:,1+ctrlshift_xy:end,:);
        original = Isyn{c1}(:,1+ctrlshift_xy:end,:);
        correlations_ctrl(c1,c2) = corr(original(:),translated(:));
        intersection_trans = min(original,translated);
        union_trans = max(original,translated);
        overlaps_ctrl(c1,c2) = sum(intersection_trans(:))./sum(union_trans(:));
        
        %2D correlation
        mip1 = max(Isyn{c1},[],3);
        mip2 = max(Isyn{c2},[],3);
        mip1_trans = max(Isyn_trans{1},[],3);
        correlations_2D(c1,c2) = corr(mip1(:),mip2(:));
        correlations_2D_ctrl(c1,c2) = corr(mip1_trans(:),mip2(:));
    end
end

%populate results structure
imdata.correlations = correlations + correlations';
imdata.overlaps = overlaps + overlaps';
imdata.correlations_ctrl = correlations_ctrl + correlations_ctrl';
imdata.overlaps_ctrl = overlaps_ctrl + overlaps_ctrl';
imdata.correlations_2D = correlations_2D + correlations_2D';
imdata.correlations_2D_ctrl = correlations_2D_ctrl + correlations_2D_ctrl';
end


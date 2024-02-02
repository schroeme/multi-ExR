function [distance_errors] = DG_ExSeq_reg_quality_vMS(img1,img2,params)
% Originally written by Dan Goodwin, modified by Margaret Schroeder
% Last update 12/14/22 (with comments), although original modifications
% made far earlier in 2022
     
    XRANGE = params.xrange;
    YRANGE = params.yrange;
    
    %determine z-range based on slices with at least some signal
    ZRANGE = [];
    for zz = 1:size(img2,3)
        slice = img2(:,:,zz);
        if ~all(slice(:) == 0)
            ZRANGE = [ZRANGE zz];
        end
    end

    img1 = img1(XRANGE,YRANGE,ZRANGE);
    img2 = img2(XRANGE,YRANGE,ZRANGE);

    if params.doplot
        figure;
        imshowpair(max(img1,[],3),max(img2,[],3));
        %subplot(1,2,1); imagesc(max(img1,[],3)); colormap gray
        %subplot(1,2,2); imagesc(max(img2,[],3)); colormap gray
    end

    %Choose N random positions throughout the volume
    N = params.N;
    xy_vol_half = params.subvol_dim/2;
    z_vol_half = min(floor(size(img1,3)/2)-1,floor(xy_vol_half*(params.xystep/params.zstep)));

    %Create the subvolumes where applicable

    %Add impossible offset so we can confirm that we're loading the
    %offsets_total correctly
    offsets_total = zeros(N,3)-30; 
    nxcorr_plot_x = zeros(N,2*xy_vol_half+1);
    nxcorr_plot_y = zeros(N,2*xy_vol_half+1);
    nxcorr_plot_z = zeros(N,2*z_vol_half+1);

    % ignores = [];

    img_offset_x = zeros(size(img1));
    img_offset_y = zeros(size(img1));
    img_offset_z = zeros(size(img1));

    %What is the minimum amount of signal in the subvolume that we would allow
    %to calculate the autocorrelation?
    thresh_img1 = prctile(img1(:),params.pct_thresh);
    thresh_img2 = prctile(img2(:),params.pct_thresh);
    i = 1;
    
    if params.doplot
        img1_bin = imbinarize(img1,thresh_img1);
        img2_bin = imbinarize(img2,thresh_img2);

        figure()
        subplot(1,2,1)
        title('Round 1')
        imagesc(max(img1_bin,[],3))
        subplot(1,2,2)
        title('Round n')
        imagesc(max(img2_bin,[],3))
    end
    xpos = randi(size(img1,1),N);
    ypos = randi(size(img1,2),N);
    zpos = randi(size(img1,3),N);

    while i<=N
        %disp(i);
    %for i = 1:N
        %Generate new -- MS: why generating new random numbers here??
        xpos(i) = randi(size(img1,2),1);
        ypos(i) = randi(size(img1,2),1);
        
        if size(img1,3) < 50  %if we have a really short stack
            zpos(i) = ceil(size(img1,3)/2);
        else
            zpos(i) = randi(size(img1,3),1);
        end

        %Confirm that the random position are within bounds
        if ~( (xpos(i)>xy_vol_half) && (xpos(i)<size(img1,1)-xy_vol_half))
    %         ignores(end+1) = i;
            continue;
        elseif ~((ypos(i)>xy_vol_half) && (ypos(i)<size(img1,2)-xy_vol_half))
    %         ignores(end+1) = i;
            continue;
        elseif ~((zpos(i)>z_vol_half) && (zpos(i)<size(img1,3)-z_vol_half))
    %         ignores(end+1) = i;
            continue
        end

        %Create the subvolumes
        subvolume1 = img1(xpos(i)-xy_vol_half:xpos(i)+xy_vol_half,...
            ypos(i)-xy_vol_half:ypos(i)+xy_vol_half,...
            zpos(i)-z_vol_half:zpos(i)+z_vol_half);
        subvolume2 = img2(xpos(i)-xy_vol_half:xpos(i)+xy_vol_half,...
            ypos(i)-xy_vol_half:ypos(i)+xy_vol_half,...
            zpos(i)-z_vol_half:zpos(i)+z_vol_half);

        subvec1 = subvolume1(:);
        subvec2 = subvolume2(:);
        
        subvec1_bin = imbinarize(subvec1,thresh_img1);
        subvec2_bin = imbinarize(subvec2,thresh_img2);
        %Confirm that the subvolumes don't have zeros in them (meaning we're
        %avoiding registration edges)
%         if any(subvec1==0) || any(subvec2==0)
%     %         ignores(end+1) = i;
%             continue
%         end
        %Set a threshold minimum amount of signal so we're not calculating
        %autocorrelation on microscope background noise
        if mean(subvec1_bin)<=0.01 || mean(subvec2_bin)<= 0.01
    %         ignores(end+1) = i;
            continue
        end

        %calculate the normalized cross correlation using code by an author on
        %the Mathworks code exchange
        [~,I_NCC,~]=template_matching(subvolume1,subvolume2);
        [val,idx] = max(I_NCC(:));
        [x,y,z] = ind2sub(size(I_NCC),idx);
        offsets_total(i,:) = [x,y,z] - ceil(size(I_NCC)/2);

        %Create 3D image volumes of all the offsets as we see it as a functino
        %of position, that way we can see if there is a 
        img_offset_x(xpos(i)-1:xpos(i)+1,ypos(i)-1:ypos(i)+1,zpos(i)-1:zpos(i)+1) = abs(offsets_total(i,1));
        img_offset_y(xpos(i)-1:xpos(i)+1,ypos(i)-1:ypos(i)+1,zpos(i)-1:zpos(i)+1) = abs(offsets_total(i,2));
        img_offset_z(xpos(i)-1:xpos(i)+1,ypos(i)-1:ypos(i)+1,zpos(i)-1:zpos(i)+1) = abs(offsets_total(i,3));


        %Note the x, y and z plots of the norm xcorrelation with the max
        %intensity projections
        test_image_mip= max(I_NCC,[],3);
        nxcorr_plot_y(i,:) = max(test_image_mip,[],1); %horizontal coord
        nxcorr_plot_x(i,:)=  max(test_image_mip,[],2); %vertical coord
        nxcorr_plot_z(i,:)= squeeze(max(max(I_NCC,[],1),[],2));%depth coord

        if mod(i,100)==0
            fprintf('Processed %i/%i\n',i,N);
        end

        %Increment the counter
        i = i+1;
    end

    %Create RGB images of the offsets, then max project them in Z 
    % img_offset_x = img_offset_x - min(img_offset_x(:));
    % img_offset_y = img_offset_y - min(img_offset_y(:));
    % img_offset_z = img_offset_z - min(img_offset_z(:));

    img_offset_x_mip = max(img_offset_x,[],3);
    img_offset_y_mip = max(img_offset_y,[],3);
    img_offset_z_mip = max(img_offset_z,[],3);

    clear rgb;
    rgb(:,:,1) = img_offset_x_mip;
    rgb(:,:,2) = img_offset_y_mip;
    rgb(:,:,3) = img_offset_z_mip;
    rgb = rgb./max(rgb(:));
    grayim = mat2gray(rgb);
    if params.doplotrgb
        figure; imshow(rgb)
        figure; imagesc(grayim)
    end
    % offsets_total(ignores,:) = [];

    N_actual = size(offsets_total,1);

    if params.doplot
        figure;
        subplot(1,3,1);
        histogram(offsets_total(:,1));
        title('Error in X (pixels)')
        subplot(1,3,2);
        histogram(offsets_total(:,2));
        title('Error in Y (pixels)')
        subplot(1,3,3);
        histogram(offsets_total(:,3));
        title('Error in Z (pixels)')
    end

    distance_errors = sqrt((params.xystep*offsets_total(:,1)).^2 + ...
        (params.xystep*offsets_total(:,2)).^2 + ...
        (params.zstep*offsets_total(:,3)).^2);
    %[Optional]Cap distance errors to it's own 99th percentile
    % distance_errors(distance_errors>prctile(distance_errors,99)) = prctile(distance_errors,99);
    if params.doplot
        figure;
        h = histogram(distance_errors,[0:.1:1]);
        h.LineWidth = 2; 
        h.FaceColor =  [0    0.4471    0.7412];
        h.FaceAlpha = 1.0;
    
        xlabel('Distance in um')
        ylabel('Count')
        title(sprintf('%s: Offsets of %i random positions in um, %0.2f within 0.5um',params.parentfolder, N_actual,sum(distance_errors<.5)/N_actual))
    end
    % saveas(gcf,sprintf('SuppFig1_MovingRun_%i.eps',moving_run),'epsc')

    %recalculate the non-normed RGB, for later comparison between the two
    %conditions
%     clear rgb;
%     rgb(:,:,1) = img_offset_x_mip;
%     rgb(:,:,2) = img_offset_y_mip;
%     rgb(:,:,3) = img_offset_z_mip;
% 
%     save(sprintf('displacement_rgb_%s.mat',params.parentfolder),'rgb')
% 
% 
%     %% Load the two displacement RGB images and then plot them in the same color space
% 
%     load(sprintf('displacement_rgb_%s.mat',params.parentfolder),'rgb');
%     rgb_gel1 = rgb;
%     filename = fullfile(imgdir,sprintf('%s_round%03d_%s.tif',params.parentfolder,fixed_run,'summedNorm'));
%     img_gel1 = load3DTif_uint16(filename);
%     XRANGE = 900:1200;
%     YRANGE = 950:1250;
%     img_gel1 = img_gel1(XRANGE,YRANGE,:);
% 
%     params.parentfolder = 'gel4N1';
%     load(sprintf('displacement_rgb_%s.mat',params.parentfolder),'rgb');
%     rgb_gel4 = rgb;
%     filename = fullfile(imgdir,sprintf('%s_round%03d_%s.tif',params.parentfolder,fixed_run,'summedNorm'));
%     img_gel4 = load3DTif_uint16(filename);
%     XRANGE = 800:1100;
%     YRANGE = 950:1250;
%     img_gel4 = img_gel4(XRANGE,YRANGE,:);
% 
%     r1 = squeeze(rgb_gel1(:,:,1)); r2 = squeeze(rgb_gel4(:,:,1));  max_r = max([r1(:); r2(:)]);
%     g1 = squeeze(rgb_gel1(:,:,2)); g2 = squeeze(rgb_gel4(:,:,2));  max_g = max([g1(:); g2(:)]);
%     b1 = squeeze(rgb_gel1(:,:,3)); b2 = squeeze(rgb_gel4(:,:,3));  max_b = max([b1(:); b2(:)]);
% 
%     rgb_gel1(:,:,1) = rgb_gel1(:,:,1)./max_r;
%     rgb_gel1(:,:,2) = rgb_gel1(:,:,2)./max_g;
%     rgb_gel1(:,:,3) = rgb_gel1(:,:,3)./max_b;
% 
%     rgb_gel4(:,:,1) = rgb_gel4(:,:,1)./max_r;
%     rgb_gel4(:,:,2) = rgb_gel4(:,:,2)./max_g;
%     rgb_gel4(:,:,3) = rgb_gel4(:,:,3)./max_b;
% 
%     fprintf('Max displacements: Y (red): %i, X (green): %i, Z (blue): %i\n',...
%         max_r, max_g, max_b);
% 
%     figure;
%     imshow(rgb_gel1); truesize; 
%     figure;
%     imshow(rgb_gel4); truesize; 
% 
%     figure;
%     subplot(1,2,1); imagesc(max(img_gel1,[],3)); 
%     subplot(1,2,2); imshow(rgb_gel1);
%     subplot(1,2,1); colormap gray; truesize; axis off;
% 
%     figure;
%     subplot(1,2,1); imagesc(max(img_gel4,[],3)); 
%     subplot(1,2,2); imshow(rgb_gel4);
%     subplot(1,2,1); colormap gray; truesize; axis off;
end

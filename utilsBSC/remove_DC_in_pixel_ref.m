function RF_PIXEL_PROCESSED = remove_DC_in_pixel_ref(MEAN_PIXEL,RF_PIXEL,WINDOW)

RF_PIXEL_PROCESSED = zeros(size(RF_PIXEL,1), size(RF_PIXEL,2), size(RF_PIXEL,3));
% RF_PIXEL_PROCESSED = zeros(size(RF_PIXEL));
O   = ones(1,size(RF_PIXEL,1));
    for i = 1 : size(RF_PIXEL,3)
        
        MAT_MEAN = O'*squeeze(MEAN_PIXEL(:,i)');
        RF_PIXEL_PROCESSED(:,:,i) = (RF_PIXEL(:,:,i) - MAT_MEAN).*WINDOW;
        clear MAT_MEAN;
    end
end
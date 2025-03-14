function [ACS_map_seg] = drawlayers(n_lay, iUS, x, z, BmodeWidth, dr)
%% ------------------------------------------------------------------------
% DESCRIPTION
% Provides:
% * ACS_map_seg     Segmented attenuation coefficient map

fig = figure;
imagesc(x, z, iUS, dr); colormap gray; axis image;
title("Do not close the figure");

ACS_map_seg = zeros(size(iUS));
layers = cell(n_lay-1,1);
for i = 1:n_lay-1
    temp_layer = drawpolyline;
%     wait(temp_layer);
    layers{i} = temp_layer;
end

prompt = "Continue? [Enter]: ";
input(prompt,"s");

for i = 1:n_lay-1
    h_layer = layers{i};
    
    h_layer.Position = sortrows(h_layer.Position,1);

    if h_layer.Position(1,1) > 0
        h_layer.Position = [0, h_layer.Position(1,2); h_layer.Position];
    end
    
    if h_layer.Position(end,1) < BmodeWidth*1e-1
        h_layer.Position = [h_layer.Position; BmodeWidth*1e-1, h_layer.Position(end,2)];
    end
    
    bw = createMask(h_layer);
    
    [~,Y] = meshgrid(1:size(bw,2),1:size(bw,1));
    
    
    [~, line_thres]=max(bw,[],1);
    
    ACS_map_seg((Y<=line_thres) & (ACS_map_seg == 0)) = i;
end
ACS_map_seg((Y>line_thres) & (ACS_map_seg == 0)) = 3;
close(fig)
end
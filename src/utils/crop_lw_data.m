function [header, data] = crop_lw_data(header, data, x_max)

% crop up to maximum time or frequency 
x_max_idx = min(round((x_max - header.xstart) / header.xstep) + 1, ...
                header.datasize(end)); 

data = data(:, :, :, :, :, 1:x_max_idx); 

header.datasize = size(data); 

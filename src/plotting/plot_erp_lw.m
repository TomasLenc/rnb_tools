function plot_erp_lw(header, data, varargin)

fs = 1/header.xstep; 

x = squeeze(data); 

if ndims(x) > 2; 
    error('only possible to plot 1D and 2D data...got %dD instead', ndims(x));  
end

plot_erp(x, fs, varargin{:}); 
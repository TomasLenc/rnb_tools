function [snippets, freq_snip_hz, freq_snip_relative] = get_fft_snips(mX, freq, frex, f_around, varargin)
% Get snippets of the spectrum around a set of center frequenices. 
% 
% Parameters
% ----------
% mX : array_like, shape=[..., frequency]
%     Raw (not noise-subtracted!!!) magnitude spectra with frequency as the 
%     last dimension. 
% freq : array_like
%     Frequencies for the FFT. 
% frex : array_like
%     Frequencies to take the snippets around. 
% f_around : float
%     The frequency range around the center frequency that will be cut out.
% 
% Returns 
% -------
% snips : array_like, shape=shape(mX)(1:end-1)
%     Average snippet around all frex. 
% idx_snip : array_like
%     Relative index of each bin in the mean_snip


parser = inputParser; 

addParameter(parser, 'avg_sides', false); 

parse(parser, varargin{:});

avg_sides = parser.Results.avg_sides; 


%%

frex_idx = dsearchn(ensure_col(freq), ensure_col(frex)); 

x_step = freq(2) - freq(1); 
n_bins_around = round(f_around / x_step); 

relative_idx = [-n_bins_around : +n_bins_around]; 

shape = size(mX); 
snippets = nan([shape(1:end-1), length(relative_idx), length(frex)]); 

for i_f=1:length(frex)
    
    idx_start = frex_idx(i_f) - n_bins_around; 
    idx_end = frex_idx(i_f) + n_bins_around; 
    
    if idx_start <= 0
        error('Frex %d is too close to the start of the input (idx_start = %d).', ...
              i_f, idx_start); 
    end
    if idx_end > shape(end)
        error('Frex %d is too close to the end of the input (idx_end = %d).', ...
              i_f, idx_end); 
    end

    index_snip = repmat({':'}, 1, ndims(mX) + 1);
    index_snip{ndims(mX) + 1} = i_f;
    
    index_mX = repmat({':'}, 1, ndims(mX));
    index_mX{end} = [idx_start : idx_end];

    snippets(index_snip{:}) = mX(index_mX{:}); 

end


if avg_sides
    
    % slice out the negative segment
    index_left = repmat({':'}, 1, ndims(mX) + 1);
    index_left{ndims(mX)} = [1 : n_bins_around+1];    
    snippets_left = snippets(index_left{:}); 

    % slice out the positive segment
    index_right = repmat({':'}, 1, ndims(mX) + 1);
    index_right{ndims(mX)} = [n_bins_around+1 : 2*n_bins_around+1];    
    snippets_right = snippets(index_right{:}); 
    
    % flip the left snippets
    snippets_left = flip(snippets_left, ndims(mX)); 

    % if the signal is complex, we need to take the complex conjugate
    % before averaging 
    if ~isreal(snippets)
        warning('the input is complex, taking complex conjugate of the left snip before averaging sides!'); 
        snippets_left = conj(snippets_left);         
    end
    
    % average the left and right snippets 
    snippets = mean(cat(ndims(snippets)+1, snippets_left, snippets_right), ...
                    ndims(snippets)+1); 
    
    freq_snip_relative = [0 : n_bins_around]; 

else
    
    freq_snip_relative = [-n_bins_around : +n_bins_around]; 

end

freq_snip_hz = freq_snip_relative * x_step; 


%% sanity check plots 

% mean_snip = mean(snippets, ndims(mX) + 1); 
% figure
% plot(freq_snip_hz, mean(mean_snip, 1), '-o')
% hold on
% plot([0,0], [0,max(mean(mean_snip, 1))], 'r')
% plot([-f_around, -f_around], [0,max(mean(mean_snip, 1))], 'r')
% plot([f_around, f_around], [0,max(mean(mean_snip, 1))], 'r')
% 
% 
% mean_snip = mean(snippets_avg_sides, ndims(mX) + 1); 
% figure
% plot(freq_snip_hz, mean(mean_snip, 1), '-o')
% hold on
% plot([0,0], [0,max(mean(mean_snip, 1))], 'r')
% plot([f_around, f_around], [0,max(mean(mean_snip, 1))], 'r')
% plot([-2/0.750 , -2/0.750], [0,max(mean(mean_snip, 1))], 'r')
% 
% 

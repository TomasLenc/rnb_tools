function [dens, x] = get_kde(data, x_start, x_end, varargin)

parser = inputParser; 

addParameter(parser, 'dt', 0.0001); 
addParameter(parser, 'sigma', 0.040); % SD
addParameter(parser, 'plot', false); 
addParameter(parser, 'amplitudes', []); 

parse(parser, varargin{:}); 

dt = parser.Results.dt; 
sigma = parser.Results.sigma; 
do_plot = parser.Results.plot; 
amplitudes = parser.Results.amplitudes; 

%%

fs = 1/dt; 

x_dur = x_end - x_start; 
x_N = round(x_dur * fs); 
x = [0 : x_N-1] / fs + x_start; 

data_idx = dsearchn(x', data); 
data_vec = zeros(1, x_N); 
if ~isempty(amplitudes)
    assert(numel(amplitudes) == numel(data_idx))
    data_vec(data_idx) = amplitudes; 
else
    data_vec(data_idx) = 1; 
end

x = ensure_row(x); 
data = ensure_col(data); 

%%

% % estiamte sigma 
% sigmas = nan(par.n_cond, 2); 
% for i_cond=1:par.n_cond
%     tmp = cellfun(@(x) x{i_cond}, {res.tap_ratios_all}, 'uni', 0); 
%     ratios = cat(1, tmp{:}); 
%     sigmas(i_cond,:) = var(ratios, [], 1); 
% end
% sigma = mean(sigmas, 1); 
% sigma = sqrt(sigma(1)); 

%%

dens = 1/(sqrt(2*pi)*sigma) * exp(-1/2 * ((data-x)/sigma).^2); 

if ~isempty(amplitudes)
    dens = dens .* amplitudes; 
end

dens = sum(dens, 1) ./ length(data); 

if do_plot
    figure
    plot(x, data_vec)
    hold on 
    plot(x, dens ./ max(dens), 'linew', 3)
    xlim([x_start, x_end])        
end

% % this is a probability distribution (but it doesn't always sum to exact
% % one especially with wider kernels, where the KDE is above zero out of range)
% assert(abs(sum(dens * dt) - 1) < 1e-6)
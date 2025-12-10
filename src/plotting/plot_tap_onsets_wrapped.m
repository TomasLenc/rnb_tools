function [f, pnl, cbar] = plot_tap_onsets_wrapped(tap_onsets_all, cycle_dur, n_cycles, varargin)
% This function plots tap onsets times as points across cycles of a
% repeating pattern.
% 
% Parameters
% ----------
% tap_onsets_all : cell of array_likes
%   One dimensional cell array. Each cell contains an array of tap onset
%   times from a particular subject or trial.
% cycle_dur : float
%   Duration of one cycle. 
% n_cycles : int
%   Total number of cycles. 
% amplitudes :  cell of array_likes, optional
%   Amplitudes corresponding to each tap onset (can be e.g. amount of force
%   estimated from a tapping sensor).
% 
% Returns 
% -------
% f : figure object 
%     Handle to the generated figure. 
% ax : axes object 
%     Handle to the generated axes. 
% cbar : colobar object 
%     Handle to the generated colobar. 
% 

parser = inputParser(); 

addParameter(parser, 'amplitudes', []); 
addParameter(parser, 'norm_force_level', 'global'); % global or local
addParameter(parser, 'stim_s', []); 
addParameter(parser, 'stim_fs', []); 

parse(parser, varargin{:}); 

amplitudes_all = parser.Results.amplitudes; 
norm_force_level = parser.Results.norm_force_level; 
stim_s = parser.Results.stim_s; 
stim_fs = parser.Results.stim_fs; 

%%

% open figure and pack subplots
f = figure('color','white',...
           'position', [67 578 343 138]); 
pnl = panel(f); 

pnl.pack('v', [25, 75]); 

ax_stim = pnl(1).select(); 
ax_tap = pnl(2).select(); 
hold(ax_stim,'on')
hold(ax_tap,'on')

cbar = []; 


if ~isempty(amplitudes_all)
    n_cols_cmap = 100; 
    cmap = brewermap(n_cols_cmap, 'Purples'); 
end

%% plot stimulus 
 
% stim_height = 0.5 * n_cycles; 
% stim_bottom = - 0.25 * n_cycles; 
stim_bottom = 0; 
stim_height = 1; 
stim_top = stim_bottom + stim_height; 

if ~isempty(stim_s) && ~isempty(stim_fs) 
    
    stim_s = stim_s - min(stim_s); 
    stim_s = stim_s ./ max(abs(stim_s)); 
    stim_s = stim_s .* stim_height + stim_bottom; 

    stim_t = [0 : length(stim_s)-1] / stim_fs; 
    
    stim_h = plot(ax_stim, stim_t, stim_s, 'linew', 0.5, 'color', [.5 .5 .5]); 
    stim_h.Color(4) = 0.7; 
end


%% plot tap onsets

for i_trial=1:length(tap_onsets_all)

    tap_onsets = tap_onsets_all{i_trial};     

    tap_onsets_mod = mod(tap_onsets, cycle_dur); 

    cycle_idx = floor(tap_onsets ./ cycle_dur); 
    
    
    if ~isempty(amplitudes_all) 
        
        % Modulate color by amplitude 
        % ---------------------------

        amplitudes = amplitudes_all{i_trial}; 
        
        if strcmp(norm_force_level, 'global')
            % normalize force by the global maximum
            amplitudes = amplitudes ./ max(1, amplitudes_all{:}); 
        elseif strcmp(norm_force_level, 'local')
            % normalize force for this trial 
            amplitudes = amplitudes ./ max(amplitudes); 
        end
                
        % prepare colors 
        edges = [0:1/n_cols_cmap:1]; 
        edges(end) = 1+1/n_cols_cmap; 
        col_idx = discretize(amplitudes, edges); 
        col = cmap(col_idx,:); 
        
        tap_onsets_h = scatter(ax_tap, tap_onsets_mod, cycle_idx, ...
                    15, 'filled'); 
        tap_onsets_h.CData = col; 
        
    else
        
        % Just tap onsets (old way)
        % -------------------------
        tap_onsets_h = scatter(ax_tap, tap_onsets_mod, cycle_idx, 15, 'filled', ...
                    'MarkerEdgeAlpha', 0.2, ...
                    'MarkerFaceAlpha', 0.2); 
        col = repmat([0.2471, 0, 0.4902], length(tap_onsets_mod), 1); 
        tap_onsets_h.CData = col; 

%         tap_onsets_h = plot(ax, tap_onsets_mod, cycle_idx, ...
%             'o','MarkerFaceColor', [0.2471    0.0000    0.4902],...
%             'MarkerEdgeColor', [0.2471    0.0000    0.4902], ...
%             'MarkerSize', 2);
%         tap_onsets_h.Color(4) = 0.2; 
    end
    
end

ax_tap.YLim = [-1, n_cycles+1]; 
ax_tap.XLim = [0, cycle_dur]; 

ax_stim.YLim = [stim_bottom, stim_top]; 
ax_stim.XLim = [0, cycle_dur]; 

ax_tap.YAxis.Visible = 'off'; 
ax_tap.XAxis.Visible = 'off'; 
ax_stim.YAxis.Visible = 'off'; 
ax_stim.XAxis.Visible = 'off'; 

pnl.de.margin = 0; 


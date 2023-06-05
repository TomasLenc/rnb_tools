function plot_erp_windows(ax, x, fs, start_times, end_times, color)

y_min = min(x); 
y_max = max(x); 
t = [0 : length(x)-1] / fs;

if length(start_times) ~= length(end_times)
    error('the number of start_times does not equal the number of end_times');
end

for i_win=1:length(start_times)

    win_start_idx = round(start_times(i_win) * fs) + 1; 
    win_end_idx = round(end_times(i_win) * fs) + 1; 

    t_win = t(win_start_idx : win_end_idx); 
    N = length(t_win); 

    fill(ax, [t_win, flip(t_win)], ...
        [repmat(y_min, 1, N), repmat(y_max, 1, N)], ...
        color, 'EdgeColor', 'none', 'FaceAlpha', 0.2); 

end



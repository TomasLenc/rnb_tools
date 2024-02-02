function tests = test_z_snr
    tests = functiontests(localfunctions);
end


function test_get_max_signif_harm_lw(test_case)
    
    n_cond = 1;
    n_trials = 1;
    trial_dur = 60;
    fs = 1000;

    [header, data] = make_lw_dataset(...
        'n_cond', n_cond, ...
        'n_trials', n_trials,...
        'trial_dur', trial_dur, ...
        'fs', fs, ...
        'data_max', 1000);
    
    [header, data] = segment_safe(header, data, {'1'}, ...
                                  'x_start', 0, 'x_duration', trial_dur);
                              
    [header, data] = RLW_FFT(header, data);
    
    n_harm_avg = 1;

    f0 = 3;
    chan_idx        = [2, 9, 11, 33];
    max_signif_harm = [1, 13, 18, 49]; 
    nonsignif_harm  = [
                        {[]} 
                        {[1]}
                        {[2,3, 10]}
                        {[20,21, 23, 30,31]}
                      ];
        
    freq_to_ignore_idx = round([f0 : f0 : fs/2] / header.xstep) + 1;
    
    data = interp_noise_bins(data, freq_to_ignore_idx, 3, 13); 
    
    for i_chan=1:length(chan_idx)
        
        signif_freq = f0 * [1 : max_signif_harm(i_chan)];
        
        signif_freq(nonsignif_harm{i_chan}) = [];
        
        signif_harm_idx = round(signif_freq / header.xstep) + 1;
        
        data(1, chan_idx(i_chan), 1, 1, 1, signif_harm_idx) = 1e9;
    end
    
%     freq = [0 : header.datasize(end)-1] * header.xstep;
%     figure
%     plot(freq, squeeze(data(1, 11, 1, 1, 1, :)))
    
    max_signif_freq_estim = get_max_signif_harm_lw(...
                            header,...
                            data, ...
                            f0, ...
                            n_harm_avg,...
                            [3, 13], ...
                            'thr_p', 0.001, ...
                            'n_successive_nonsignif_steps_to_stop', 3);      
                                
    % check that we found all the channels
    assert(all(find(max_signif_freq_estim) == chan_idx'));
    
    % check that the max signif freq for chosen channels was estimated
    % correctly
    max_signif_freq_true = f0 * max_signif_harm;

    assert(all(max_signif_freq_estim(max_signif_freq_estim > 0) == ...
               max_signif_freq_true'));
    
end


function test_get_max_signif_harm_lw_avg3(test_case)
    
    n_cond = 1;
    n_trials = 1;
    trial_dur = 60;
    fs = 1000;

    [header, data] = make_lw_dataset(...
        'n_cond', n_cond, ...
        'n_trials', n_trials,...
        'trial_dur', trial_dur, ...
        'fs', fs, ...
        'data_max', 1000);
    
    [header, data] = segment_safe(header, data, {'1'}, ...
                                  'x_start', 0, 'x_duration', trial_dur);
                              
    [header, data] = RLW_FFT(header, data);
    
    n_harm_avg = 3;
    
    f0 = 3;
    chan_idx        = [2, 9, 11, 33];
    max_signif_harm = [1, 13, 25, 49]; 
    nonsignif_harm  = [
                        {[]} 
                        {[1]}
                        {[2,3,4, 10]}
                        {[20,21, 23, 30,31]}
                      ];
        
    freq_to_ignore_idx = round([f0 : f0 : fs/2] / header.xstep) + 1;
    
    data = interp_noise_bins(data, freq_to_ignore_idx, 3, 13); 
    
    for i_chan=1:length(chan_idx)
        
        signif_freq = f0 * [1 : max_signif_harm(i_chan)];
        
        signif_freq(nonsignif_harm{i_chan}) = [];
        
        signif_harm_idx = round(signif_freq / header.xstep) + 1;
        
        data(1, chan_idx(i_chan), 1, 1, 1, signif_harm_idx) = 1e9;
    end
    
%     freq = [0 : header.datasize(end)-1] * header.xstep;
%     figure
%     plot(freq, squeeze(data(1, 11, 1, 1, 1, :)))
    
    max_signif_freq_estim = get_max_signif_harm_lw(...
                            header,...
                            data, ...
                            f0, ...
                            n_harm_avg,...
                            [3, 13], ...
                            'thr_p', 0.001, ...
                            'n_successive_nonsignif_steps_to_stop', 2);      
                                
    % check that we found all the channels
    assert(all(find(max_signif_freq_estim) == chan_idx'));
    
    % check that the max signif freq for chosen channels was estimated
    % correctly
    max_signif_freq_true = f0 * (floor(max_signif_harm / n_harm_avg) * n_harm_avg + n_harm_avg);

    assert(all(max_signif_freq_estim(max_signif_freq_estim > 0) == ...
               max_signif_freq_true'));
    
end










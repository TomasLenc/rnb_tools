function tests = test_phase_locking
    tests = functiontests(localfunctions);
end


function test_get_phase_locking_noise(test_case)
    
    trial_dur = 600;
    fs = 500;

    [header, data] = make_dataset(...
        'n_cond', 1, ...
        'n_chans', 2, ...
        'n_trials', 3,...
        'trial_dur', trial_dur, ...
        'fs', fs, ...
        'fill_with', 'noise');
    
    [header, data] = RLW_segmentation(header, data, {'1'}, ...
        'x_start', 0, 'x_duration', trial_dur);
    
    [r, theta, ph] = get_phase_locking(data, 1.25, fs, 'skip_s_buffer', 20);
    
%     figure
%     polarplot(ph(:), ones(numel(ph), 1), 'o',...
%               'color',[89, 89, 89]/255)
%     hold on
%     polarplot(repmat(theta(1), 1, 2), [0, r(1)],...
%               'color', 'red', ...
%               'linew',1.5)

    % check that the resulting vector length is small 
    assert(all(r(:) < 0.1));

end


function test_get_phase_locking_perfect_periodic(test_case)
    
    fs = 500;
    
    ir = get_square_kernel(fs, ...
        'duration', 0.100, ...
        'rampon', 0, ...
        'rampoff', 0 ...
        ); 
    
    % make clean signal for the whole trial 
    [x, t] = get_s(...
                    [1 0 0 0 1 0 0 0 1 0 0 0], ...
                    0.2, ...
                    fs, ...
                    'n_cycles', 32, ...
                    'ir', ir ...
                    );
                
    [r, theta, ph] = get_phase_locking(x, 1.25, fs, 'skip_s_buffer', 10);
    
    % check that the resulting vector is equal to 1
    assert(abs(1 - r) < 1e-10);

    
end


function test_get_phase_locking_multiD(test_case)
    
    fs = 500;
    
    ir = get_square_kernel(fs, ...
        'duration', 0.100, ...
        'rampon', 0, ...
        'rampoff', 0 ...
        ); 
    
    % make clean signal for the whole trial 
    [x_trial, t] = get_s(...
                    [1 0 0 0 1 0 0 0 1 0 0 0], ...
                    0.2, ...
                    fs, ...
                    'n_cycles', 32, ...
                    'ir', ir ...
                    );
                
    x = []; 
    x(:, 1, 1, 1, 1, :) = repmat(x_trial, [10, 1]); 
    
    [r, theta, ph] = get_phase_locking(x, 1.25, fs, 'skip_s_buffer', 10);
    
    % check that the resulting vector is equal to 1
    assert(length(unique(r)) == 1);
    assert(length(unique(theta)) == 1);

    
end












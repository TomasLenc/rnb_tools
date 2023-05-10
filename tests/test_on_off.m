function tests = test_on_off
    tests = functiontests(localfunctions);
end


function test_get_on_off_zeros(test_case)
    
    trial_dur = 60;
    fs = 500;

    x = zeros(10, 64, 1, 1, 1, round(trial_dur * fs)); 
    
    [on_off_contrast] = get_on_off(x, fs, 0.8, 0.050);
    
    % check that the resulting vector length is small 
    assert(all(isnan(on_off_contrast(:))));

end


function test_get_on_off_perfect_periodic(test_case)
    
    fs = 1000;
    
    ir = get_square_kernel(fs, ...
        'duration', 0.025, ...
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
                
    [on_off_contrast, on_val, off_val] = get_on_off(x, fs, 0.8, 0.050);
    
    assert(on_val == 0.5);
    assert(off_val == 0);
    assert(on_off_contrast == 1);
    
end


function test_get_on_off_perfect_equal(test_case)
    
    fs = 500;
    
    ir = get_square_kernel(fs, ...
        'duration', 0.050, ...
        'rampon', 0, ...
        'rampoff', 0 ...
        ); 
    
    % make clean signal for the whole trial 
    [x, t] = get_s(...
                    [1 0 1 0 1 0 1 0 1 0 1 0], ...
                    0.2, ...
                    fs, ...
                    'n_cycles', 32, ...
                    'ir', ir ...
                    );
                
    [on_off_contrast] = get_on_off(x, fs, 0.8, 0.050);
    
    assert(on_off_contrast == 0);
    
end


function test_get_on_off_perfect_offbeat(test_case)
    
    fs = 500;
    
    ir = get_square_kernel(fs, ...
        'duration', 0.050, ...
        'rampon', 0, ...
        'rampoff', 0 ...
        ); 
    
    % make clean signal for the whole trial 
    [x, t] = get_s(...
                    [0 0 1 0 0 0 1 0 0 0 1 0], ...
                    0.2, ...
                    fs, ...
                    'n_cycles', 32, ...
                    'ir', ir ...
                    );
                
    [on_off_contrast] = get_on_off(x, fs, 0.8, 0.050);
    
    assert(on_off_contrast == -1);
    
end


function test_get_on_off_emphasis(test_case)
    
    fs = 500;
    
    ir = get_square_kernel(fs, ...
        'duration', 0.100, ...
        'rampon', 0, ...
        'rampoff', 0 ...
        ); 
    
    emph_magns = [0, 0.1, 0.5, 1]; 
    
    res = nan(1, length(emph_magns)); 
    
    for i_emph=1:length(emph_magns)
        
        % make clean signal for the whole trial 
        [x_trial, t] = get_s(...
                        [1 1 1 1 0 1 1 1 0 0 1 0], ...
                        0.2, ...
                        fs, ...
                        'n_cycles', 32, ...
                        'ir', ir, ...
                        'emph_period', 4, ...
                        'emph_magn', emph_magns(i_emph) ...
                        );

        x = []; 
        x(:, 1, 1, 1, 1, :) = repmat(x_trial, [10, 1]); 

        [on_off_contrast] = get_on_off(x, fs, 0.8, 0.050);
    
        assert(length(unique(on_off_contrast)) == 1);
        
        res(i_emph) = on_off_contrast(1); 
        
    end
    
    % we should get increasing values of the on/off contrast with increasing
    % on-beat emphasis
    assert(all(sort(res) == res))
    
end












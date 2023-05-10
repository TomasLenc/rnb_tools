function tests = test_phase_locking
    tests = functiontests(localfunctions);
end


%     figure
%     polarplot(ph(:), ones(numel(ph), 1), 'o',...
%               'color',[89, 89, 89]/255)
%     hold on
%     polarplot(repmat(theta(1), 1, 2), [0, r(1)],...
%               'color', 'red', ...
%               'linew',1.5)


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


function test_get_phase_locking_perfect_periodic_multiD(test_case)
    
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
    
    % check that all the results are equal to 1
    assert(length(unique(r)) == 1);
    assert(length(unique(theta)) == 1);

    
end


function test_get_phase_locking_jitter(test_case)
    
    fs = 500;
    
    ir = get_square_kernel(fs, ...
        'duration', 0.100, ...
        'rampon', 0, ...
        'rampoff', 0 ...
        ); 
    
    jitters = [0, 0.020, 0.050, 0.070]; 
    
    r = nan(1, length(jitters));
    
    for i_jit=1:length(jitters)
        
        % make signal for the whole trial 
        [x, t] = get_s(...
                        [1 0 0 0 1 0 0 0 1 0 0 0], ...
                        0.2, ...
                        fs, ...
                        'n_cycles', 32, ...
                        'ir', ir, ...
                        'jitter', jitters(i_jit) ...
                        );

        r(i_jit) = get_phase_locking(x, 1.25, fs, 'skip_s_buffer', 10);
        
    end
    
    % check that r is gradually decreasing with increaasing amount of jitter
    assert(all(sort(r, 'descend') == r))
    
end









function tests = test_circular
    tests = functiontests(localfunctions);
end


function test_phase_from_onset_times_clean(test_case)

    n_taps = 100;
    period = 0.6;
    noise_sd = period * 0.0;
    onset_times = [0 : n_taps-1] * period + randn(1, n_taps) * noise_sd;

    plausible_periods = 0.6;

    [ph, r, mean_ph, closest_period] = phase_from_onset_times(...
                                    onset_times, plausible_periods...
                                    );

    assert(all(ph == 0));
    assert(r == 1);
    assert(mean_ph == 0);
    assert(closest_period == period);

end
    

function test_phase_from_onset_times_find_period(test_case)

    n_taps = 100;
    period = 0.6;
    noise_sd = period * 0.1;
    
    plausible_periods = [0.4, 0.6, 0.8];

    while 1
        onset_times = [0 : n_taps-1] * period + randn(1, n_taps) * noise_sd;
        [~, idx] = min(abs(median(diff(onset_times)) - plausible_periods));
        if idx == 2
            break
        end
    end

    [ph, r, mean_ph, closest_period] = phase_from_onset_times(...
                                    onset_times, plausible_periods...
                                    );

    assert(closest_period == period);

end
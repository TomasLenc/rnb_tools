function tests = test_segmentation
    tests = functiontests(localfunctions);
end


function test_epoch_chunks_1D(test_case)

    fs = 44100; 
    n_cycles = 16; 
    ir = get_erp_kernel(fs); 
    x = get_s([1,1,1,0,1,1,1,0,1,1,0,0], 0.2, fs, ...
              'n_cycles', n_cycles, 'ir', ir); 
    x_chunked = epoch_chunks(x, fs, 2.4); 
    
    assert(all(size(x_chunked) == [n_cycles, round(fs * 0.2*12)])); 
    assert(all(x_chunked(1, :) == x_chunked(end, :))); 
   
end


function test_epoch_chunks_noise_fail(test_case)

    fs = 44100; 
    n_cycles = 16; 
    ir = get_erp_kernel(fs); 
    x = randn(1, round(n_cycles * 12 * 0.2 * fs)); 
    x_chunked = epoch_chunks(x, fs, 2.4); 
    
    assert(all(size(x_chunked) == [n_cycles, round(fs * 0.2*12)])); 
    assert(~all(x_chunked(1, :) == x_chunked(end, :))); 
   
end


function test_epoch_chunks_1D_offset_start(test_case)

    fs = 44100; 
    n_cycles = 16; 
    ir = get_erp_kernel(fs); 
    x = get_s([1,1,1,0,1,1,1,0,1,1,0,0], 0.2, fs, ...
              'n_cycles', n_cycles, 'ir', ir); 
    x_chunked = epoch_chunks(x, fs, 2.4, 'start', 2.4 * 2); 
    
    assert(all(size(x_chunked) == [n_cycles - 2, round(fs * 0.2*12)])); 
    assert(all(x_chunked(1, :) == x_chunked(end, :))); 
   
end


function test_epoch_chunks_multiD(test_case)

    fs = 44100; 
    n_cycles = 16; 
    
    ir = get_erp_kernel(fs); 
    x1 = get_s([1,1,1,0,1,1,1,0,1,1,0,0], 0.2, fs, ...
              'n_cycles', n_cycles, 'ir', ir); 
          
    ir = get_square_kernel(fs); 
    x2 = get_s([1,1,1,0,1,1,1,0,1,1,0,0], 0.2, fs, ...
              'n_cycles', n_cycles, 'ir', ir); 
          
    x = [];
    x(1, 1, :, :) = [x1; x1]; 
    x(1, 2, :, :) = [x2; x2]; 
  
    x_chunked = epoch_chunks(x, fs, 2.4); 
    
    assert(all(size(x_chunked) == [n_cycles, 1, 2, 2, round(fs * 0.2*12)]))
    
    assert(all(x_chunked(1, 1, 1, 1, :) == x_chunked(end, 1, 1, 1, :))); 
    assert(all(x_chunked(1, 1, 2, 1, :) == x_chunked(end, 1, 2, 1, :))); 
    assert(~all(x_chunked(1, 1, 1, 1, :) == x_chunked(1, 1, 2, 1, :))); 
   
end


function test_segment_safe_one_cond(test_case)
    % Test epoching with one trigger code value. 
    
    n_cond = 3;
    n_trials = 5;
    trial_dur = 10;
    pause_dur = 3;

    [header, data] = make_lw_dataset(...
        'n_cond', n_cond, ...
        'n_trials', n_trials,...
        'trial_dur', trial_dur, ...
        'pause_dur', pause_dur,...
        'fs', 1000);

    [header_ep, data_ep] = segment_safe(...
        header, data, {'1'}, 'x_start', 0, 'x_duration', trial_dur);
    
    assert(isequal([header_ep.events.epoch], [1, 2, 3, 4, 5]));
    assert(all([header_ep.events.latency] == 0));
    assert(all(strcmp({header_ep.events.code}, '1')));
    
    assert(all(data_ep(:, 1, 1, 1, 1, 1) == 1));

end


function test_segment_safe_multi_cond(test_case)
    % Test epoching with multiple trigger code values. 
    n_cond = 3;
    n_trials = 5;
    trial_dur = 10;
    pause_dur = 3;

    [header, data] = make_lw_dataset(...
        'n_cond', n_cond, ...
        'n_trials', n_trials,...
        'trial_dur', trial_dur, ...
        'pause_dur', pause_dur,...
        'fs', 1000);

    [header_ep, data_ep] = segment_safe(...
        header, data, {'1', '2'}, 'x_start', 0, 'x_duration', trial_dur);
    
    assert(isequal([header_ep.events.epoch], [1 : 2 * n_trials]));
    assert(all([header_ep.events.latency] == 0));
    assert(all(cellfun(@(x) str2num(x), {header_ep.events.code}) == ...
               data_ep(:, 1, 1, 1, 1, 1)'));
end


function test_segment_safe_multi_cond_overlap(test_case)
    % Test epoching with multiple trigger code values, and the epoch duration
    % is so long that other trigger values end up in the epoch. The function
    % should delete any additional triggers that have been used for segmentation
    %  and keep the rest. 
    n_cond = 3;
    n_trials = 5;
    trial_dur = 10;
    pause_dur = 3;

    [header, data] = make_lw_dataset(...
        'n_cond', n_cond, ...
        'n_trials', n_trials,...
        'trial_dur', trial_dur, ...
        'pause_dur', pause_dur,...
        'end_recording_buffer', 10, ...
        'fs', 1000);

    [header_ep, data_ep] = segment_safe(...
        header, data, {'1', '2'},...
        'x_start', 0, ...
        'x_duration', trial_dur + pause_dur + 1);
    
    % check that latency of events used for segmentation is always 0 (all the
    % other overlapping evets have been deleted)
    assert(all(ismember(...
        {header_ep.events([header_ep.events.latency] == 0).code}, {'1', '2'} ...
        )));
    % check that the additional events are kept in the header
    assert(all(strcmp(...
        {header_ep.events([header_ep.events.latency] ~= 0).code}, '3' ...
        )));
        
%     % visual check 
%     struct2table(header.events)
%     struct2table(header_ep.events)
    
end


function test_segment_safe_too_long(test_case)
    % Test the case when the requested epoch goes outside of data range. Should
    % raise an exception. 
    n_cond = 1;
    n_trials = 5;
    trial_dur = 10;
    pause_dur = 3;

    [header, data] = make_lw_dataset(...
        'n_cond', n_cond, ...
        'n_trials', n_trials,...
        'trial_dur', trial_dur, ...
        'pause_dur', pause_dur,...
        'end_recording_buffer', 0.010, ...
        'fs', 1000);

    verifyError(test_case, ...
        @()  segment_safe(...
            header, data, {'1'}, ...
            'x_start', 0, ...
            'x_duration', trial_dur + 0.5), ...
        'segment_safe:SegmentOutOfDataRange');
    
end


function test_segment_safe_too_long_ignore(test_case)
    % Test the case when the requested epoch goes outside of data range. Should
    % raise an exception. 
    n_cond = 1;
    n_trials = 5;
    trial_dur = 10;
    pause_dur = 3;

    [header, data] = make_lw_dataset(...
        'n_cond', n_cond, ...
        'n_trials', n_trials,...
        'trial_dur', trial_dur, ...
        'pause_dur', pause_dur,...
        'end_recording_buffer', 0.010, ...
        'fs', 1000);

    [header_ep, data_ep] = segment_safe(...
        header, data, {'1'}, ...
        'x_start', 0, ...
        'x_duration', trial_dur + 0.5, ...
        'out_of_range', 'ignore');
    
    assert(header_ep.datasize(1) == 4);
    assert(isequal([header_ep.events.epoch], [1 : n_trials-1]));
    assert(all([header_ep.events.latency] == 0));
    assert(all(cellfun(@(x) str2num(x), {header_ep.events.code}) == ...
               data_ep(:, 1, 1, 1, 1, 1)'));
        
end


function test_segment_safe_too_long_zero_pad(test_case)
    % Test the case when the requested epoch goes outside of data range (it
    % ends after the data ends). This can happen e.g. when we ask for too
    % much buffer after trial end but the acquisition has been stopped. We
    % will zero pad the missing data here. 
    n_cond = 1;
    n_trials = 5;
    trial_dur = 10;
    pause_dur = 3;

    [header, data] = make_lw_dataset(...
        'n_cond', n_cond, ...
        'n_trials', n_trials,...
        'trial_dur', trial_dur, ...
        'pause_dur', pause_dur,...
        'end_recording_buffer', 0.010, ...
        'fs', 1000);

    [header_ep, data_ep] = segment_safe(...
        header, data, {'1'}, ...
        'x_start', 0, ...
        'x_duration', trial_dur + 0.5, ...
        'out_of_range', 'zero_pad');
    
    assert(header_ep.datasize(1) == n_trials);
    assert(isequal([header_ep.events.epoch], [1 : n_trials]));
    assert(all([header_ep.events.latency] == 0));

    % make sure the data value at time 0 is exactly equal to the event code
    % in that epoch 
    assert(all(cellfun(@(x) str2num(x), {header_ep.events.code}) == ...
               data_ep(:, 1, 1, 1, 1, 1)'));
     
    % find indices where the padding should be       
    expected_N_pad = (header.events(end).latency / header.xstep + (trial_dur + 0.5) / header.xstep) - ...
                      header.datasize(end); 
                  
   % ensure the padding is in the last epoch 
   assert(all(data_ep(end, 1, 1, 1, 1, end-expected_N_pad+1 : end) == 0))
   
   % in the last epoch, the value just before padding shouldn't be zero 
   assert(data_ep(end, 1, 1, 1, 1, end-expected_N_pad) ~= 0)

   % all other epochs shouldn't have any padding
   for i=1:n_trials-1
        assert(~all(data_ep(i, 1, 1, 1, 1, end-expected_N_pad+1 : end) == 0))
   end
   
 
end


function test_segment_safe_too_late_zero_pad(test_case)
    % Test the case when the requested epoch starts before the continuous
    % data started. This can happen e.g. when we ask for too much buffer
    % before trial onset. We will zero pad the missing data here.
    n_cond = 1;
    n_trials = 5;
    trial_dur = 10;
    pause_dur = 3;

    [header, data] = make_lw_dataset(...
        'n_cond', n_cond, ...
        'n_trials', n_trials,...
        'trial_dur', trial_dur, ...
        'pause_dur', pause_dur,...
        'start_recording_buffer', 0.010, ...
        'end_recording_buffer', 0.010, ...
        'fs', 1000);

    [header_ep, data_ep] = segment_safe(...
        header, data, {'1'}, ...
        'x_start', -0.5, ...
        'x_duration', trial_dur + 0.5, ...
        'out_of_range', 'zero_pad');
    
    assert(header_ep.datasize(1) == n_trials);
    assert(isequal([header_ep.events.epoch], [1 : n_trials]));
    assert(all([header_ep.events.latency] == 0));
    
    % make sure the data value at time 0 is exactly equal to the event code
    % in that epoch 
    idx_time_zero = 0.5 / header.xstep + 1; 
    assert(all(cellfun(@(x) str2num(x), {header_ep.events.code}) == ...
               data_ep(:, 1, 1, 1, 1, idx_time_zero)'));
     
    % find indices where the padding should be       
    expected_N_pad = (0.5 - 0.010) / header.xstep; 
    
    % ensure the padding is in the first epoch 
    assert(all(data_ep(1, 1, 1, 1, 1, 1:expected_N_pad) == 0))
                 
    % in the first epoch, the value just after padding shouldn't be zero 
    assert(data_ep(1, 1, 1, 1, 1, expected_N_pad+1) ~= 0)
    
    % all other epochs shouldn't have any padding
    for i=2:n_trials
        assert(~all(data_ep(i, 1, 1, 1, 1, 1:expected_N_pad) == 0))
    end

 
end


function test_segment_safe_negative_buffer(test_case)
    % Test when x_start is a negative value. 
    n_cond = 1;
    n_trials = 5;
    trial_dur = 10;
    pause_dur = 6;
    start_recording_buffer = 8;
    end_recording_buffer = 8;

    [header, data] = make_lw_dataset(...
        'n_cond', n_cond, ...
        'n_trials', n_trials,...
        'trial_dur', trial_dur, ...
        'pause_dur', pause_dur,...
        'start_recording_buffer', start_recording_buffer, ...
        'end_recording_buffer', end_recording_buffer, ...
        'fs', 1000);
    
    ep_buffer = 5;

    [header_ep, data_ep] = segment_safe(...
        header, data, {'1'},...
        'x_start', -ep_buffer, ...
        'x_duration', trial_dur + 2 * ep_buffer);
    
    assert(isequal([header_ep.events.epoch], [1 : n_trials]));
    assert(all(strcmp({header_ep.events.code}, '1')));
    assert(all([header_ep.events.latency] == 0));
    assert(header_ep.xstart == -ep_buffer);
    assert(header_ep.datasize(6) * header_ep.xstep == trial_dur + 2 * ep_buffer);
    
    idx0 = round(ep_buffer / header_ep.xstep + 1); 
    assert(all(data_ep(:, 1, 1, 1, 1, idx0) == 1));

end


function test_segment_safe_resegment(test_case)
    
    n_cond = 3;
    n_trials = 5;
    trial_dur = 10;
    pause_dur = 6;
    start_recording_buffer = 8;
    end_recording_buffer = 8;

    [header, data] = make_lw_dataset(...
        'n_cond', n_cond, ...
        'n_trials', n_trials,...
        'trial_dur', trial_dur, ...
        'pause_dur', pause_dur,...
        'start_recording_buffer', start_recording_buffer, ...
        'end_recording_buffer', end_recording_buffer, ...
        'fs', 1000);
    
    ep_buffer = 5;

    [header_ep, data_ep] = segment_safe(...
        header, data, {'1', '3'},...
        'x_start', -ep_buffer, ...
        'x_duration', trial_dur + 2 * ep_buffer);
    
    [header_ep, data_ep] = segment_safe(...
        header_ep, data_ep, {'1', '3'},...
        'x_start', 0, ...
        'x_duration', trial_dur);
    
    assert(isequal([header_ep.events.epoch], [1 : 2 * n_trials]));
    assert(all([header_ep.events.latency] == 0));
    assert(all(cellfun(@(x) str2num(x), {header_ep.events.code}) == ...
               data_ep(:, 1, 1, 1, 1, 1)'));
end



function test_segment_safe_resegment_too_long(test_case)
    
    n_cond = 3;
    n_trials = 5;
    trial_dur = 10;
    pause_dur = 6;
    start_recording_buffer = 8;
    end_recording_buffer = 8;

    [header, data] = make_lw_dataset(...
        'n_cond', n_cond, ...
        'n_trials', n_trials,...
        'trial_dur', trial_dur, ...
        'pause_dur', pause_dur,...
        'start_recording_buffer', start_recording_buffer, ...
        'end_recording_buffer', end_recording_buffer, ...
        'fs', 1000);
    
    [header_ep, data_ep] = segment_safe(...
        header, data, {'1', '3'},...
        'x_start', 0, ...
        'x_duration', trial_dur);
    
    verifyError(test_case, ...
        @()  segment_safe(...
            header_ep, data_ep, {'1', '3'},...
            'x_start', 0, ...
            'x_duration', trial_dur + 3), ...
        'segment_safe:SegmentOutOfDataRange');

end


function test_segment_safe_resegment_too_long_ignore(test_case)

    trial_dur = 10;
    pause_dur = 2;
    start_recording_buffer = 5;
    end_recording_buffer = 0;
    ep_buffer = 1;

    [header, data] = make_lw_dataset(...
        'event_sequence', {'1', '2', '3', '4'}, ...
        'trial_dur', trial_dur, ...
        'pause_dur', pause_dur,...
        'start_recording_buffer', start_recording_buffer, ...
        'end_recording_buffer', end_recording_buffer, ...
        'fs', 1000);
    
    [header_ep, data_ep] = segment_safe(...
        header, data, {'1', '3', '4'},...
        'x_start', -ep_buffer, ...
        'x_duration', trial_dur + pause_dur + 2, ...
        'out_of_range', 'ignore');
    
    % we should have 3 events in the output
    assert(length(header_ep.events) == 3);
    
    % but only two epochs in the data
    assert(size(data_ep, 1) == 2);
    
    % first, there should be event code "1". Everything fine. 
    assert(strcmp(header_ep.events(1).code, '1'))
    assert(header_ep.events(1).latency == 0)
    assert(header_ep.events(1).epoch == 1)

    % But for the first epoch, there was an overlaping event with code "2".
    % Now we didn't ask to segment based on code "2" so the event will be
    % registered in the output
    assert(strcmp(header_ep.events(2).code, '2'))
    assert(header_ep.events(2).latency == trial_dur + pause_dur)
    assert(header_ep.events(2).epoch == 1)
        
    % The second epoch is based on the following trigger "3" and even
    % though there is an overlap with the following trigger "4" within the
    % duration of the epoch, that one is ignored since it's one of the
    % codes used for the current segmentation (the function "knows" this is
    % just the following epoch)
    assert(strcmp(header_ep.events(3).code, '3'))
    assert(header_ep.events(3).latency == 0)
    assert(header_ep.events(3).epoch == 2)
    
    % another check - the simualtino function placed an impulse in the data
    % at the start of each trial - we can make sure we got the correct data
    % samples to start the trial 
    idx0 = round(ep_buffer / header_ep.xstep + 1); 
    assert(isequal(data_ep(:, 1, 1, 1, 1, idx0), [1; 3]));

    % re-epoch, this time include all event codes, but ignore out-of-range
    % epochs -> this should remove the event '2' that was in the first
    % epoch 
    [header_ep, data_ep] = segment_safe(...
        header_ep, data_ep, {'1', '2', '3', '4'},...
        'x_start', 0, ...
        'x_duration', trial_dur-6, ...
        'out_of_range', 'ignore');
    
    assert(strcmp(header_ep.events(1).code, '1'))
    assert(header_ep.events(1).latency == 0)
    assert(header_ep.events(1).epoch == 1)
    
    assert(strcmp(header_ep.events(2).code, '3'))
    assert(header_ep.events(2).latency == 0)
    assert(header_ep.events(2).epoch == 2)
       
    assert(length(header_ep.events) == 2);
    assert(size(data_ep, 1) == 2);
    
    assert(isequal(data_ep(:, 1, 1, 1, 1, 1), [1; 3]));
        
end



function test_segment_safe_additional_event_fields(test_case)
    % Same test rationale as above, but this time we create an additional
    % field in the events structure of the intput header. The field should
    % be carreid into the output header after epoching. 
    
    trial_dur = 10;
    pause_dur = 2;
    start_recording_buffer = 5;
    end_recording_buffer = 0;
    ep_buffer = 1;

    [header, data] = make_lw_dataset(...
        'event_sequence', {'1', '2', '3', '4'}, ...
        'trial_dur', trial_dur, ...
        'pause_dur', pause_dur,...
        'start_recording_buffer', start_recording_buffer, ...
        'end_recording_buffer', end_recording_buffer, ...
        'fs', 1000);
    
    header.events(1).polarity = 'a';  
    header.events(2).polarity = 'b';  
    header.events(3).polarity = 'c';  
    header.events(4).polarity = 'd';  
    
    [header_ep, data_ep] = segment_safe(...
        header, data, {'1', '3', '4'},...
        'x_start', -ep_buffer, ...
        'x_duration', trial_dur + pause_dur + 2, ...
        'out_of_range', 'ignore');
    
    assert(isequal(...
        header_ep.events(1),...
        struct('code', '1', 'latency', 0, 'epoch', 1, 'polarity', 'a', ...
               'latency_orig', start_recording_buffer)...
        )); 
    
    assert(isequal(...
        header_ep.events(2), ...
        struct('code', '2', 'latency', trial_dur + pause_dur, 'epoch', 1, 'polarity', 'b', ...
               'latency_orig', start_recording_buffer + trial_dur + pause_dur)...
        ));     
    
    assert(isequal(...
        header_ep.events(3), ...
        struct('code', '3', 'latency', 0, 'epoch', 2, 'polarity', 'c', ...
               'latency_orig', start_recording_buffer + 2*(trial_dur + pause_dur))...
        )); 
    
end
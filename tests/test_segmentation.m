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

    [header, data] = make_dataset(...
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

    [header, data] = make_dataset(...
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

    [header, data] = make_dataset(...
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
end


function test_segment_safe_too_long(test_case)
    % Test the case when the requested epoch goes outside of data range. Should
    % raise an exception. 
    n_cond = 1;
    n_trials = 5;
    trial_dur = 10;
    pause_dur = 3;

    [header, data] = make_dataset(...
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

    [header, data] = make_dataset(...
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
        'ignore_out_of_range', true);
    
    assert(header_ep.datasize(1) == 4);
    assert(isequal([header_ep.events.epoch], [1 : n_trials-1]));
    assert(all([header_ep.events.latency] == 0));
    assert(all(cellfun(@(x) str2num(x), {header_ep.events.code}) == ...
               data_ep(:, 1, 1, 1, 1, 1)'));
        
end


function test_segment_safe_negative_buffer(test_case)
    % Test when x_start is a negative value. 
    n_cond = 1;
    n_trials = 5;
    trial_dur = 10;
    pause_dur = 6;
    start_recording_buffer = 8;
    end_recording_buffer = 8;

    [header, data] = make_dataset(...
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

    [header, data] = make_dataset(...
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

    [header, data] = make_dataset(...
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

    [header, data] = make_dataset(...
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
        'ignore_out_of_range', true);
    
    assert(isequal(...
        header_ep.events(1),...
        struct('code', '1', 'latency', 0, 'epoch', 1)...
        )); 
    assert(isequal(...
        header_ep.events(2), ...
        struct('code', '2', 'latency', trial_dur + pause_dur, 'epoch', 1)...
        ));     
    assert(isequal(...
        header_ep.events(3), ...
        struct('code', '3', 'latency', 0, 'epoch', 2)...
        )); 
    
    assert(length(header_ep.events) == 3);
    assert(size(data_ep, 1) == 2);
    
    idx0 = round(ep_buffer / header_ep.xstep + 1); 
    assert(isequal(data_ep(:, 1, 1, 1, 1, idx0), [1; 3]));

    % re-epoch
    [header_ep, data_ep] = segment_safe(...
        header_ep, data_ep, {'1', '2', '3', '4'},...
        'x_start', 0, ...
        'x_duration', trial_dur-6, ...
        'ignore_out_of_range', true);
    
    assert(isequal(...
        header_ep.events(1),...
        struct('code', '1', 'latency', 0, 'epoch', 1)...
        )); 
    assert(isequal(...
        header_ep.events(2), ...
        struct('code', '3', 'latency', 0, 'epoch', 2)...
        )); 
    
    assert(length(header_ep.events) == 2);
    assert(size(data_ep, 1) == 2);
    
    assert(isequal(data_ep(:, 1, 1, 1, 1, 1), [1; 3]));
        
end


% 
% 
% importLetswave(6);
% if ~isdir('tmp')
%     mkdir('tmp');
% end
% header.name = 'test_dataset';
% CLW_save('tmp', header, data);
% % CLW_save('tmp', header_ep, data_ep);
% importLetswave(7);
% cd tmp
% letswave
% 
% cd ..
% rmdir('tmp', 's');
% importLetswave(6);
% 
% 

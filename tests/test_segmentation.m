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

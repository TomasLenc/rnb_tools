function tests = test_noise_gen
    tests = functiontests(localfunctions);
end


function test_col_noise_gen(test_case)

    fs = 1000; 

    dur = 200; 

    N = fs * dur; 

    noise_exp = -1.5; 

    x = get_colored_noise(N, fs, noise_exp); 

    freq = [0 : N/2] / N * fs; 

    [pxx, f] = pwelch(x, fs); 

    f = log(f(5:end)); 
    pxx = log(pxx(5:end)); 

    c = polyfit(f, pxx, 1);

%     plot(f, pxx); 
%     c(1)

    assert(abs(c(1) - noise_exp) < 0.1)

end










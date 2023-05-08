function [r, theta, ph] = get_phase_locking(x, f, fs, varargin)
% Calculate phase stability of narrowband-filtered signal phase across pulse 
% positions. 
% 
% Parameters
% ----------
% x : array_like, shape=[..., time]
%     Time-domain signal to analyse. Time must be the last dimension. 
% f : float
%     Center frequency for the narrowband band-pass filter. 
% fs : float
%     Sampling rate (samples/s)
% sigma : float, optional, default=f/2
%     Standard deviation of the gaussian frequency domain response of the kernel
% pulse_period_sec : float, optional, default=1/f
%     Period of the pulse that is going to be analysed (in seconds). 
% pulse_phase_sec : float, optional, default=0
%     Phase of the pulse that is going to be analysed (in seconds). This is with
%     respect to the start of the signal (i.e. time of the first pulse position 
%     after begining of the signal).   
% data_type : str, optional, {'continuous', 'discrete'}, default='continuous'
%     Type of data that is passed in. 
% plot_ir : bool, optional, default=false
%     If true, a diagnostic figure of the impulse response kernel is generated.
% skip_s_buffer : float, optional, default=0
%     Number of seconds at the begining and the end of the signal that will
%     likely contain filtering artifacts. These segments will be discarded. 
% 
% Returns 
% -------
% r : float
%     Length of the average phase vector. 
% theta : float
%     Angle of the average phase vector. 

parser = inputParser; 

addParameter(parser, 'sigma', f/2); 
addParameter(parser, 'plot_ir', false); 
addParameter(parser, 'data_type', 'continuous'); 
addParameter(parser, 'pulse_period_sec', 1/f); 
addParameter(parser, 'pulse_phase_sec', 0); 
addParameter(parser, 'skip_s_buffer', 0); 

parse(parser, varargin{:});

sigma = parser.Results.sigma; 
plot_ir = parser.Results.plot_ir; 
data_type = parser.Results.data_type; 
pulse_period_sec = parser.Results.pulse_period_sec; 
pulse_phase_sec = parser.Results.pulse_phase_sec; 
skip_s_buffer = parser.Results.skip_s_buffer; 


%% phase-locking analysis

% band-pass filter the data
cmw = get_wavelet_kernel(f, sigma, fs, 'do_plot', plot_ir);

% Depending on whether the data is discrete or continuous, we have two
% approaches to calculating the phase. 
if strcmpi(data_type, 'discrete')
%     % In case of discrete process, we can just take the onset-times of the 
%     % impulses and consider those points to represent 2pi phase of the system. 
%     % This is equivalent to working with asynchronies which is typical in
%     % the SMS literature (see Repp). 
%     discrOnsets = discrOnsetsAllCond{iCond}; 
%     quantTapOnsets = round((discrOnsets-pulsePhaseTime)/pulsePeriodTime)*pulsePeriodTime + pulsePhaseTime; 
%     asy = discrOnsets - quantTapOnsets; 
%     sdAsy =  std(asy); 
%     % convert asynchronies to circular (i.e. rescale to 0-2pi)
%     ph = asy/pulsePeriodTime*2*pi; 
    
elseif strcmpi(data_type, 'continuous')
    
    n_cmw = length(cmw); 
    hn_cmw = floor(n_cmw / 2) + 1;
    n_data = size(x, ndims(x)); 
    n_conv = n_data + n_cmw - 1; 

    % fft of wavelet
    cmwX = fft(cmw, n_conv); 
    % normalize amplitude
    cmwX = cmwX ./ max(cmwX); 
    % fft of data
    dataX = fft(x, n_conv, ndims(x)); 
    
    % convolution
    % -----------   
    
    % time is on the last dimension - we will loop over everything else
    
    % allocate 
    shape = size(x); 
    convers = nan([shape(1:end-1), n_conv]); 
    
    % prepare index vector as cell, e.g. {1,1,1,1,':'} 
    nv = ndims(x) - 1;  % exclude last dimension
    ready = false; 
    idx_while_loop = [repmat({1}, 1, nv), {':'}]; 
    while ~ready
      
        convers(idx_while_loop{:}) = ifft(...
                       ensure_row(cmwX) .* ...
                       ensure_row(squeeze(dataX(idx_while_loop{:})))...
                       ); 
                        
        % Update the index vector:
        % Assume that the WHILE loop is ready
        ready = true;       
         % Loop through dimensions
        for k = 1:nv        
            % Increase current index by 1
            idx_while_loop{k} = idx_while_loop{k} + 1;   
            % Does it exceed the length of current dim?
            if idx_while_loop{k} <= shape(k) 
               % No, WHILE loop is not ready now
               ready = false;  
               % idx_while_loop(k) increased successfully, leave "for k" loop
               break;         
            end
            % Reset idx_while_loop{k}, proceed to next k
            idx_while_loop{k} = 1;          
        end
    end    
  
    index = repmat({':'}, 1, ndims(x)); 
    index{end} = [hn_cmw : size(convers, ndims(x))-hn_cmw+1]; 
    convers = convers(index{:});
    
    % In case of continuous process, we have to estimate the phase
    % (which is a mess if we have multicomponent impulse response!). 
    % A typical way to do this is filter-hilbert. 
    ph = angle(convers); 
    
    % Then we want to see what the phase of the system was at the
    % metric pulse positions (we would like to see consistent phase value 
    % if the system was synchronized). 
    x_dur = n_data / fs; 
    pulse_pos_sec = [0 + pulse_phase_sec : pulse_period_sec : x_dur-1/fs]; 
    
    % remove start and end buffer that may contain filtering artifacts
    pulse_pos_sec(pulse_pos_sec < skip_s_buffer) = []; 
    pulse_pos_sec(pulse_pos_sec > x_dur-skip_s_buffer) = []; 
        
    pulse_pos_idx = round(pulse_pos_sec * fs) + 1; 
    
    index = repmat({':'}, 1, ndims(x)); 
    index{end} = pulse_pos_idx; 
    ph = ph(index{:}); 
    
else
    
    error('data type "%s" not supported', data_type); 
    
end

% get mean vector length and mean phase
index = repmat({':'}, 1, ndims(x)); 
index{end} = 1; 
index(shape == 1) = {1};

r = []; 
r(index{:}) = abs(mean(exp(1j * ph), ndims(x))); 

theta = []; 
theta(index{:}) = angle(mean(exp(1j * ph), ndims(x))); 
















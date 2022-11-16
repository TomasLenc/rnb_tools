function ir = get_square_kernel(fs, varargin)
% Generate square-wave impulse response kernel. 
% 
% Parameters
% ----------
% fs : int
%     Sampling rate in Hz. 
% duration : float, default=0.1
%     Duration of the up state in seconds. 
% rampon : float, default=0
%     Duration of the linear onset ramp in seconds. 
% rampoff : float, default=0
%     Duration of the linear offset ramp in seconds. 
% 
% Returns
% -------
% ir : array_like
%     Time-domain kernel. 
%


parser = inputParser; 

addParameter(parser, 'duration', 0.1, @isfloat); 
addParameter(parser, 'rampon', 0, @isfloat); 
addParameter(parser, 'rampoff', 0, @isfloat); 

parse(parser, varargin{:}); 

duration = parser.Results.duration; 
rampon = parser.Results.rampon; 
rampoff = parser.Results.rampoff; 


ir = ones(1, round(duration * fs));

% apply onset and offset ramp
ir(1 : round(rampon * fs)) = ir(1 : round(rampon * fs)) .* ...
                            linspace(0, 1, round(rampon * fs));

ir(end - round(rampoff * fs) + 1 : end) = ir(end - round(rampoff * fs) + 1 : end) .* ...
                                     linspace(1, 0, round(rampoff * fs));


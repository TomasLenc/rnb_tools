function [ir, t_ir] = get_erp_kernel(fs, varargin)
% Generate erp-like impulse response kernel. 
% 
% Parameters
% ----------
% fs : int
%     Sampling rate in Hz. 
% duration : float, default=0.2
%     Duration of the whole kernel. 
% f0s : array of floats, default=7.
%     Frequencies of the constituting sinusoidal components in Hz. 
% t0s : float, default=0
%     Time offsets of each constituting sinusiodal component in s. 
% taus : array of floats, default=0.050
%     Decay constants of the constituting sinusiodal components in s. 
% amplitudes : array of floats, default=1
%     Amplitudes of the constituting sinusiodal components.  
% 
% Returns
% -------
% ir : array_like
%     Time-domain kernel. 
%

parser = inputParser; 

addParameter(parser, 'amplitudes', 1); 
addParameter(parser, 'duration', 0.2, @isfloat); 
addParameter(parser, 'f0s', 7); 
addParameter(parser, 't0s', 0); 
addParameter(parser, 'taus', 0.050); 

parse(parser, varargin{:}); 

duration = parser.Results.duration; 
f0s = parser.Results.f0s; 
t0s = parser.Results.t0s; 
taus = parser.Results.taus; 
amplitudes = parser.Results.amplitudes; 


N_ir = round(duration * fs); 
t_ir = [0 : N_ir - 1] / fs; 

ir = zeros(1, N_ir);
for fi=1:length(f0s)

    ir = ir + amplitudes(fi) * ...
             (t_ir - t0s(fi)) / taus(fi) .* ...
             exp( 1 - (t_ir - t0s(fi)) / taus(fi) ) .* ...
             sin( 2 * pi * f0s(fi) * (t_ir - t0s(fi)) ) ; 

end
ir = ir ./ length(f0s); 

% normalize to integreate to 1
ir = ir ./ sum(ir); 


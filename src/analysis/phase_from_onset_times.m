function [ph, r, mean_ph, closest_period] = phase_from_onset_times(...
                            onset_times, plausible_periods, varargin)
% This function takes onset times of a point process timeseries (e.g. tapping, 
% where each time point is a defined as a contact between the finger and the 
% tapping surface). It calculates phase angles with respect to a requested 
% target period. It returns the phase angles and also a mean angle and mean
% vector length. Phase 0 of the target period is assumed to be at time 0.  
% 
% Parameters
% ----------
% onset_times : array_like, 1D
%     One-dimensional vector of onset times. 
% plausible_periods : array_like
%     One or more period(s) to which the tapping (or other analysed point 
%     process) could be synchronized. Must in in the same units as onset_times. 
%     If only a single value is passed, it will be used directly to calculate
%     phase. If multiple values are passed, the one closest to the median
%     inter-onset interval is chosen. 
% 
% Returns 
% -------
% ph : array_like
%     Vector with phase angle assigned to each onset time. 
% r : float
%     Mean vector length. 
% mean_ph : float
%     Phase angle of the mean vector. 
% closest_period : float
%     If multiple plausible periods were passed, this is the one closest to the 
%     median inter-onset interval (i.e. phases were calculated wrt this period).
% 

% get a vector of inter-onset intervals
iois = diff(sort(onset_times));
% pick plausible period that is the closest to the median inter-onset interval
[~, closest_period_idx] = min(abs(median(iois) - plausible_periods));
closest_period = plausible_periods(closest_period_idx);

% get phase angles
ph = mod(onset_times, closest_period) / closest_period * 2 * pi;

mean_vector = mean(exp(1i*ph));
mean_ph = angle(mean_vector);
r = abs(mean_vector);



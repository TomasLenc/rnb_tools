function [r, y_pred, y_test_out, x_test_out] = ...
            run_trf_encoding_loc(x_all_trials, y_all_trials, fs, par_trf)
% leave-one-trial out crossvalidation - pass x and y as cell arrays 

assert(size(y_all_trials,1) == size(x_all_trials,1), ...
    'y has %d samples and X has %d', size(y_all_trials,1), size(x_all_trials,1))

n_cv_folds = length(x_all_trials); 

r = nan(1, n_cv_folds); 

for i_test_fold=1:n_cv_folds

    test_idx = i_test_fold; 
    
    train_idx = ones(1, n_cv_folds); 
    train_idx(test_idx) = 0; 
    train_idx = find(train_idx); 
    
    x_train = x_all_trials(train_idx); 
    y_train = y_all_trials(train_idx); 
    
    x_test = x_all_trials(test_idx); 
    y_test = y_all_trials(test_idx); 
    
    % cross validate labda
    cv = mTRFcrossval(x_train, y_train, fs, ...
                      par_trf.direction, ...
                      par_trf.tmin, ...
                      par_trf.tmax, ...
                      par_trf.lambdas, ...
                      'fast', 0,...
                      'method', 'Tikhonov', ...
                      'verbose', 0);

    % Get optimal hyperparameters
    [rmax, idx] = max(mean(mean(cv.r), 3));
    lambda = par_trf.lambdas(idx);

    model = mTRFtrain(x_train, y_train, ...
                      fs, par_trf.direction, ...
                      par_trf.tmin,...
                      par_trf.tmax, ...
                      lambda, ...
                      'method','Tikhonov', ...
                      'verbose', 0);

   [y_pred, test] = mTRFpredict(x_test, y_test, model, ...
                              'zeropad',1, 'verbose', 0);

   r(i_test_fold) = test.r; 

end


y_test_out = y_test{:}; 
x_test_out = x_test{:}; 
function clusters_clean = cluster_perm_test1d(data, varargin)
% Perform a cluster-based permutation test. First, it will do a
% one sample, sign-permutation test to identify clusters of data points
% above the requested statistical threshold. Then it will generate a null
% distribution of cluster sizes (again using the sign-permuted data).
% Finally, it will return observed clusters that are bigger than expected
% by chance (i.e. greater than what's expected from the null distribution
% of cluster sizes. 
% 
% Parameters
% ----------
% data : array_like, shape=[observations, time]
%     Data to test. It must be 2-dimensional, with observations (e.g.
%     subjects) on the first dimension. The second dimension can be
%     anything, typically time. 
% n_perm : int, optional, default=10000
%     Number of permutations
% p_thr : float, optional, default=0.05
%     Statistical threshold, pvalue, between 0 and 1. 
% p_clust : float, optional, default=0.05
%     Statistical threshold, pvalue, to decide which clusters will be
%     considered significant based on the null distribution of cluster
%     sizes. 
% mu : float, optional, default=0
%     Value (mean) to test the data against. 
% do_plot: bool, optional, default=false
%     If true, the results will be plotted. 
%     
% Returns
% -------
% clusters_clean : struct
%     Strucure with indices of each significant cluster in the field named
%     "PixelIdxList". 
% 
parser = inputParser(); 

addParameter(parser, 'n_perm', 10000); 
addParameter(parser, 'p_thr', 0.05); 
addParameter(parser, 'p_clust', 0.05); 
addParameter(parser, 'tail', 0); % -1 below, 1 above, 0 both
addParameter(parser, 'mu', 0); 
addParameter(parser, 'plot', false); 

parse(parser, varargin{:}); 

n_perm = parser.Results.n_perm; 
p_thr = parser.Results.p_thr; 
p_clust = parser.Results.p_clust; 
tail = parser.Results.tail; 
mu = parser.Results.mu; 
do_plot = parser.Results.plot; 

%% geneate null distribution 

n_sub = size(data, 1); 

% initialize
data_fake = zeros(size(data));
data_perm = zeros(n_perm, size(data, 2));

for i_perm = 1:n_perm
    % random sequence of +/-1
    rand_signs = sign(randn(n_sub, 1));
    % randomize sign of r across subjects
    for i_sub=1:n_sub
        data_fake(i_sub,:) = rand_signs(i_sub) * (data(i_sub,:) - mu);
    end
    % average across subjects for each time point 
    data_perm(i_perm,:) = mean(data_fake, 1);
end

if tail==0
    p_thr = p_thr/2; 
end
    
ci = prctile(data_perm, [100*p_thr, 100-100*p_thr], 1); 

%% test observed data against the null distribution 

data_mean_obs = mean(data - mu, 1); 

if tail==1
    signif = data_mean_obs > ci(2,:); 
elseif tail==0
    signif = data_mean_obs < ci(1,:) | data_mean_obs > ci(2,:); 
elseif tail==-1
    signif = data_mean_obs < ci(1,:); 
end

clusters = bwconncomp(logical(signif));

%% test permuted data against the null distribution to get null distr of cluster sizes

% initialize cluster sizes from permutation
fake_clust_sizes = zeros(n_perm,1);

for i_perm=1:n_perm    
    if tail==1
        signif_fake = data_perm(i_perm,:) > ci(2,:); 
    elseif tail==0
        signif_fake = data_perm(i_perm,:) < ci(1,:) | data_perm(i_perm,:) > ci(2,:); 
    elseif tail==-1
        signif_fake = data_perm(i_perm,:) < ci(1,:); 
    end
    % identify clusters
    clusters_fake = bwconncomp(logical(signif_fake));
    % find cluster sizes
    clust_ns = cellfun(@length, clusters_fake.PixelIdxList);
    if ~isempty(clust_ns)
        fake_clust_sizes(i_perm) = max(clust_ns);
    end
end

% compute cluster threshold
clust_thr = prctile(fake_clust_sizes, 100-p_clust*100);

%% test cluster sizes against null 

% keep only clusters above threshold
clusters_clean = clusters; 
clusters_clean.PixelIdxList = []; 
for i_clust=1:clusters.NumObjects
    if numel(clusters.PixelIdxList{i_clust}) > clust_thr
        clusters_clean.PixelIdxList{end+1} = clusters.PixelIdxList{i_clust}; 
    end
end
clusters_clean.NumObjects = numel(clusters.PixelIdxList); 

%% plot 

if do_plot
%     % show distribution of cluster sizes
%     figure
%     hold on
%     histogram(fake_clust_sizes)
%     plot([1 1]*clust_thr, get(gca,'ylim'),'r--','linew',3)

    % plot data 
    f = figure('color', 'white', 'pos', [156 237 542 295]); 
    pnl = panel(f); 
    ax = pnl.select(); 
    hold on 
    t = [1 : size(data, 2)]; 
    plot([1, size(data, 2)], [mu, mu], 'k', 'linew', 3); 
    plot(t, data_mean_obs+mu, 'linew', 2)
    h = fill([t, flip(t)], [mu+data_mean_obs+ci(1,:), flip(mu+data_mean_obs + ci(2,:))], 'k'); 
    h.FaceAlpha = 0.2; 
    h.EdgeColor = 'none'; 
    for i_clust=1:length(clusters_clean.PixelIdxList)
        plot(t(clusters_clean.PixelIdxList{i_clust}), ...
             repmat(mu+min(data_mean_obs)+min(ci(1,:)), size(clusters_clean.PixelIdxList{i_clust})), ...
             'r', 'linew', 3)
    end
    pnl.fontsize = 16; 
    pnl.margin = [15, 15, 10, 10]; 
end


function stat_struct = permutation_test_time2(data1,data2, cfg)

% intitialize stat sturct
stat_struct = [];
[h,p,confi,stats] = ttest(data1,data2);
stat_struct.t_wave = squeeze(stats.tstat);

%pull out parameters
time_vec = cfg.time_vec;
calc_lims = cfg.calc_lims;
alpha = cfg.alpha;
nperm = cfg.nperm;

%% calculate observed values
zero_ind = find((time_vec <= calc_lims(1)) | (time_vec >= calc_lims(2)));
alphaup = 1-cfg.alpha/2;
t_crit = tinv(alphaup,nanmean(stats.df));
tstat_vec = stats.tstat;
t_above_ind = tstat_vec >= t_crit;
CC_pos = bwconncomp(t_above_ind);
pos_comp = CC_pos.PixelIdxList;
pos_comp = pos_comp(find(cell2mat(cellfun(@(x) length(x) > 20, pos_comp,'UniformOutput',0))));
%negative
t_below_ind = tstat_vec <= -1.*t_crit;
CC_neg = bwconncomp(t_below_ind);
neg_comp = CC_neg.PixelIdxList;
neg_comp = neg_comp(find(cell2mat(cellfun(@(x) length(x) > 20, neg_comp,'UniformOutput',0))));
ind_obs = cat(2,pos_comp,neg_comp);

if isempty(ind_obs)
    stat_struct.PixelIdxList = {};
    stat_struct.t_obs = 0;
    stat_struct.p_obs = 1;
    stat_struct.PixelIdxList_all = {};
    stat_struct.t_obs_all = 0;
    stat_struct.p_obs_all = 1;
    stat_struct.specReturn = [];
    return
end


%set up structure
stat_struct.PixelIdxList = ind_obs;
stat_struct.t_obs = zeros(length(ind_obs),1);
stat_struct.p_obs = ones(length(ind_obs),1);
%get tsum
for obs = 1:length(ind_obs)
    stat_struct.t_obs(obs) = sum(tstat_vec(ind_obs{obs}));
end

maxt_perm = zeros(1,nperm);
mint_perm = zeros(1,nperm);

%set up data
n_obs = size(data1,1);

%run permutation
parfor r = 1:nperm
    rand_vals = rand(n_obs,1);
    low_ind = find(rand_vals<0.5);
    high_ind = find(rand_vals>=0.5);
    
    %stack random data
    rand_dat1 = cat(1,data1(low_ind,:),...
        data2(high_ind,:));
    rand_dat2 = cat(1,data2(low_ind,:),...
        data1(high_ind,:));
    %compute test
    [h,p,ci,stats] = ttest2(rand_dat1,rand_dat2,'Dim',1);
    tstat_vec_rand = stats.tstat;
    t_above_ind = tstat_vec_rand >= t_crit;
    CC_pos = bwconncomp(t_above_ind);
    pos_comp = CC_pos.PixelIdxList;
    pos_comp = pos_comp(find(cell2mat(cellfun(@(x) length(x) > 20, pos_comp,'UniformOutput',0))));
    if ~isempty(pos_comp)
        %get tsum
        pos_t = zeros(length(pos_comp),1);
        for obs = 1:length(pos_t)
            pos_t(obs) = sum(tstat_vec_rand(pos_comp{obs}));
        end
        maxt_perm(r) = max(pos_t);
    end
    %negative
    t_below_ind = tstat_vec_rand <= -1.*t_crit;
    CC_neg = bwconncomp(t_below_ind);
    neg_comp = CC_neg.PixelIdxList;
    neg_comp = neg_comp(find(cell2mat(cellfun(@(x) length(x) > 20, neg_comp,'UniformOutput',0))));
    if ~isempty(neg_comp)
        %get tsum
        neg_t = zeros(length(neg_comp),1);
        for obs = 1:length(neg_t)
            neg_t(obs) = sum(tstat_vec_rand(neg_comp{obs}));
        end
        mint_perm(r) = min(neg_t);
    end
end%rand perm


valid_sig = zeros(length(ind_obs),1);
for c = 1:length(ind_obs)
    if stat_struct.t_obs(c) > 0
        p_obs = sum(maxt_perm > stat_struct.t_obs(c))/nperm;
    else
        p_obs = sum(mint_perm < stat_struct.t_obs(c))/nperm;
    end
    stat_struct.p_obs(c) = p_obs;
    if p_obs < (alpha)
        valid_sig(c) = 1;
    end
end

%save all p-values for fdr correction
stat_struct.PixelIdxList_all = stat_struct.PixelIdxList;
stat_struct.t_obs_all = stat_struct.t_obs;
stat_struct.p_obs_all = stat_struct.p_obs;

valid_sig = find(valid_sig);

if isempty(valid_sig)
    stat_struct.specReturn = [];
end

stat_struct.PixelIdxList = stat_struct.PixelIdxList(valid_sig);
stat_struct.t_obs = stat_struct.t_obs(valid_sig);
stat_struct.p_obs = stat_struct.p_obs(valid_sig);

%get bounding box
specReturn = zeros(size(stat_struct.t_wave));
for c = 1:length(valid_sig)
    %significance mat
    specReturn(stat_struct.PixelIdxList{c}) = 1;
end
stat_struct.specReturn = specReturn;

end
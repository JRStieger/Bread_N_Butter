function [specReturn,sig_struct] = CalcPermSpec2(ersp1,ersp2,time_vec,calc_lims,freqs,cfg)

%first spec
%valid_height = find(~isnan(squeeze(ersp1(:,10,100))));
%ersp1 = ersp1(valid_height,:,:);
zero_ind = find((time_vec <= calc_lims(1)) | (time_vec >= calc_lims(2)));
ersp_tmp1 = ersp1;
ersp_tmp1(:,:,zero_ind) = 0;
%second_spec
%valid_height = find(~isnan(squeeze(ersp2(:,10,100))));
%ersp2 = ersp2(valid_height,:,:);
zero_ind = find((time_vec <= calc_lims(1)) | (time_vec >= calc_lims(2)));
ersp_tmp2 = ersp2;
ersp_tmp2(:,:,zero_ind) = 0;

%get observed values
[h,p,confi,stats] = ttest2(ersp_tmp1,ersp_tmp2);
alpha = cfg.alpha;
sig_spots = squeeze(p<alpha);
se = strel('disk',3,6);
sig_spots = imopen(sig_spots,se);
CC = bwconncomp(sig_spots);
if isfield(cfg,'thresh')
    thresh = cfg.thresh;
else
    thresh = 0.005*size(ersp1,2)*size(ersp1,3);
end
nperm = cfg.nperm;
con_spots = find(cell2mat(cellfun(@(x) length(x), CC.PixelIdxList,'UniformOutput',0))>thresh);

if isempty(con_spots)
    specReturn = [];
    sig_struct = [];
    return
end

sig_struct.PixelIdxList = CC.PixelIdxList(con_spots);
sig_struct.t_obs = zeros(length(con_spots),1);
sig_struct.p_obs = zeros(length(con_spots),1);

t_obs_spec = squeeze(stats.tstat);
for c = 1:length(con_spots)
    sig_struct.t_obs(c) = sum(t_obs_spec(CC.PixelIdxList{con_spots(c)}));
end

maxt_perm = zeros(1,nperm);
mint_perm = zeros(1,nperm);

dat_tot = [ersp_tmp1;ersp_tmp2];
len1 = height(ersp_tmp1);
len2 = height(ersp_tmp2);

%tic
parfor r = 1:nperm
    rand_ind = randperm(height(dat_tot));
    rand_dat1 = dat_tot(rand_ind(1:len1),:,:);
    rand_dat2 = dat_tot(rand_ind((len1+1):end),:,:);

    %compute t stats
    [h,p_rand,confi,stats_rand] = ttest2(rand_dat1,rand_dat2);
    sig_spots_rand = squeeze(p_rand<alpha);
    
    sig_spots_rand = imopen(sig_spots_rand,se);
    CC_rand = bwconncomp(sig_spots_rand);
    con_spots_rand = find(cell2mat(cellfun(@(x) length(x), CC_rand.PixelIdxList,'UniformOutput',0))>thresh);
    
    if isempty(con_spots_rand)
        continue
    end
    
    t_rand_spec = squeeze(stats_rand.tstat);
    t_rand = zeros(length(con_spots_rand),1);
    for c = 1:length(con_spots_rand)
        t_rand(c) = sum(t_rand_spec(CC_rand.PixelIdxList{con_spots_rand(c)}));
    end
    
    maxt_perm(r) = max(t_rand);
    mint_perm(r) = min(t_rand);
    
end%for permutation
%toc
valid_sig = zeros(length(con_spots),1);
for c = 1:length(con_spots)
    if sig_struct.t_obs(c) > 0
        p_obs = sum(maxt_perm > sig_struct.t_obs(c))/nperm;
    else
        p_obs = sum(mint_perm < sig_struct.t_obs(c))/nperm;
    end
    sig_struct.p_obs(c) = p_obs;
    if p_obs < (alpha/cfg.ncorrect)
        valid_sig(c) = 1;
    end
    
end

valid_sig = find(valid_sig);

if isempty(valid_sig)
    specReturn = [];
end

sig_struct.PixelIdxList = sig_struct.PixelIdxList(valid_sig);
sig_struct.t_obs = sig_struct.t_obs(valid_sig);
sig_struct.p_obs = sig_struct.p_obs(valid_sig);

%get bounding box
time_mat = repmat(time_vec,length(freqs),1);
freq_mat = repmat(freqs',1,length(time_vec));
freq_mat_ind = repmat((1:length(freqs))',1,length(time_vec));
time_lims = zeros(length(valid_sig),2);
freq_lims = zeros(length(valid_sig),2);
rec_lims = zeros(length(valid_sig),4);

specReturn = zeros(size(sig_spots));
for c= 1:length(valid_sig)
    %freq_lims
    tmp_freq = freq_mat(sig_struct.PixelIdxList{c});
    freq_lims(c,:) = [min(tmp_freq),max(tmp_freq)];
    %time_lims
    tmp_time = time_mat(sig_struct.PixelIdxList{c});
    time_lims(c,:) = [min(tmp_time),max(tmp_time)];
    %bounding box
    %freq_ind_lims
    tmp_freq = freq_mat_ind(sig_struct.PixelIdxList{c});
    tmp_lims = [min(tmp_freq),max(tmp_freq)];
    rec_lims(c,:) = [time_lims(c,1),...
        tmp_lims(1),...
        time_lims(c,2)-time_lims(c,1),...
        tmp_lims(2)-tmp_lims(1)];
    %significance mat
    specReturn(sig_struct.PixelIdxList{c}) = 1;
    
end

sig_struct.freq_lims = freq_lims;
sig_struct.time_lims = time_lims;
sig_struct.rectangle = rec_lims;

sig_struct.specReturn = specReturn;

end
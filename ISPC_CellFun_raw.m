function [ISPC] = ISPC_CellFun_raw(chan1_dat,chan2_dat)
phase_diff = chan1_dat - chan2_dat;
phase_exp = exp(1i.*(phase_diff));
ISPC = abs(nanmean(phase_exp,1));
end


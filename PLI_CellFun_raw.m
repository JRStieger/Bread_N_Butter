function [PLI] = PLI_CellFun_raw(chan1_dat,chan2_dat)
phase_diff = chan1_dat - chan2_dat;
phase_exp = sign(phase_diff);
PLI = abs(mean(phase_exp,1,'omitnan'));
end


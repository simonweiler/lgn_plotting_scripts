function [peak_red peak_blue red_irrad blue_irrad] = seq_photo(trace,photodiode,redpeak_start,redpeak_end,bluepeak_start,bluepeak_end)

BL_AMPA_P1 = [redpeak_start-28 redpeak_start-1];

BL_AMPA_P2 = [bluepeak_start-0.5 bluepeak_start-1];
RW_AMPA_P1 = [redpeak_start+3 redpeak_start+30];
RW_AMPA_P2 = [bluepeak_start+3 bluepeak_start+30];

end
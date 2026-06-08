ensureRotpendHardwarePath(mfilename('fullpath'));

cfg = rotpend_hwinit(0.001);
h = cfg.h;
daoutoffs = cfg.daoutoffs;
daoutgain = cfg.daoutgain;
adinoffs = cfg.adinoffs;
adingain = cfg.adingain;

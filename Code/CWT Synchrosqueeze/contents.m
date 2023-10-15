% ASST Toolbox
% Version 1.0 January 2019
% Copyright (c) 2017-2019 by Lin Li (China) and Qingtang Jiang (USA)
% 
%
% MAIN - This is an example for the two-component LFM signal
% RENYI_ENTROPY - Find the time-varying sigma by Renyi entropy
% CWT_MC_TI - Computes the CWT Ws for a fixed time ti
% BOUNDARY_CWT - Calculate the boundaries of multicomponent on the CWT plane;
% IMAGESQ - Displays time-frequency result of the Synchrosqueezing transform
% CWT_SST_MC_TI - Computes the CWT Ws, the phase transform of CWT Rs, and the  synchrosqueezing transform SST Ts of the signal s.
% CWT_SST_MC - Computes the CWT Ws, the phase transform of CWT Rs, and the  synchrosqueezing transform SST Ts of the signal s.
% RENYI_ENTROPY_GLOBAL - Computes concentration of a time-frquency distribution.
% CWT_2ND_ORDER_SST_TV - Computes the ADAPTIVE Continuous Wavelet transform: Ws, the ADAPTIVE 1st-order time-varying synchrosqueezing transform: Ts1_vt.
% the ADAPTIVE 2nd-order time-varying synchrosqueezing transform: Ts2_vt.   
% and the phase transform w2nd.


% This program can not be used for commercialization without the authorization of its author(s). 
% If you use this toolbox for research, please must cite the following paper:
% [1] Lin Li, Haiyan Cai and Qingtang Jiang, "Adaptive synchrosqueezing transform with a time-varying parameter for non-stationary signal separation," preprint, 2018, arXiv:1812.11364.
% [2] Lin Li, Haiyan Cai, Hongxia Han, Qingtang Jiang, and Hongbing Ji, "Adaptive short-time Fourier transform and synchrosqueezing transform for non-stationary signal separation," preprint, 2018, arXiv:1812.11292.
% [3] Haiyan Cai, Qingtang Jiang, Lin Li, and Bruce W. Suter, "Analysis of adaptive short-time Fourier transform-based synchrosqueezing", preprint, 2018, arXiv:1812.11033.
% For any comment or bug report, please send e-mail to lilin@xidian.edu.cn or jiangq@umsl.edu. 

function Re = Renyi_entropy_global(Vs,al)
% Renyi_ENTROPY_GLOBAL computes concentration of a time-frquency
% distribution.
%
%
% input  
%      Vs:   a kind of time-frquency distribution
%       al:   a constant, al>2
%
% output 
%       Re:  the  Renyi entropy of a time-frquency distribution
%
% By Dr. Lin Li, October 2017
% Copyright (c) 2017 Xidian University
% 
% 
% This program can not be used for commercialization without the authorization of its author(s). 

% If you use this toolbox for research, please must cite the following paper:
% [1] Lin Li, Haiyan Cai and Qingtang Jiang, "Adaptive synchrosqueezing transform with a time-varying parameter for non-stationary signal separation," preprint, 2018, arXiv:1812.11364.
% [2] Lin Li, Haiyan Cai, Hongxia Han, Qingtang Jiang, and Hongbing Ji, "Adaptive short-time Fourier transform and synchrosqueezing transform for non-stationary signal separation," preprint, 2018, arXiv:1812.11292.
% [3] Haiyan Cai, Qingtang Jiang, Lin Li, and Bruce W. Suter, "Analysis of adaptive short-time Fourier transform-based synchrosqueezing", preprint, 2018, arXiv:1812.11033.

% For any comment or bug report, please send e-mail to lilin@xidian.edu.cn or jiangq@umsl.edu. 



Ds = Vs;
alfa=2*al;
rb1 = sum(sum((abs(Ds(:,:))).^alfa));
rb2 = (sum(sum(abs(Ds(:,:)).^2))).^al;
Re = 1/(1-al)*log2(rb1/rb2);

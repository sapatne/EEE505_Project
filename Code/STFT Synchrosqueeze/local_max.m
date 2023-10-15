function a_1 = local_max(F1,gamma)
% LOCAL_MAX finds the local maxima
%
%  Input:
%   F1: a function
%   gamma: a threshold
%
% Output:
%   a_1: the positions of the local maxima
%
% Example:
% a_1 = local_max(F1,0.3);
%
%
% Lin Li, March 2018
% Copyright (c) 2018 Xidian University
% 
% 
% This program can not be used for commercialization without the authorization of its author(s). 

% If you use this toolbox for research, please must cite the following papers:
% [1] Lin Li, Haiyan Cai and Qingtang Jiang, "Adaptive synchrosqueezing transform with a time-varying parameter for non-stationary signal separation," preprint, 2018, arXiv:1812.11364.
% [2] Lin Li, Haiyan Cai, Hongxia Han, Qingtang Jiang, and Hongbing Ji, "Adaptive short-time Fourier transform and synchrosqueezing transform for non-stationary signal separation," preprint, 2018, arXiv:1812.11292.
% [3] Haiyan Cai, Qingtang Jiang, Lin Li, and Bruce W. Suter, "Analysis of adaptive short-time Fourier transform-based synchrosqueezing", preprint, 2018, arXiv:1812.11033.

% For any comment or bug report, please send e-mail to lilin@xidian.edu.cn
% or jiangq@umsl.edu. 


tep = abs(F1);
tep = tep(:);
tep = tep';
tep1 = tep/max(tep);
l_tep1 = [0,0,0,tep1(4:end-1),0];
r_tep1 = [0,0,0,tep1(2:end-3),0];
m_tep1 = [0,0,0,tep1(3:end-2),0];
a_1 = find(((m_tep1>l_tep1)+(m_tep1>r_tep1))==2 & m_tep1>gamma*max(abs(tep1)))-1;
end
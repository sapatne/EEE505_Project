function F1 = stft_ti(s,ci,ti)
% STFT_TI computes the STFT for fixed sigma and t 
%
%  Input:
%   s: signal
%   ci: sigma
%   ti: time instant t
%
% Output:
%   F1: STFT
%   fj: frequency bins
%
% Example:
% F1 = stft_ti(s,0.02,128);
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


ol = length(s);  % original length
s = s(:);      % Turn into column vector
s = s - mean(s);      %zero mean
n = length(s);     % checking length of signal
p2 = ceil(log2(n));
s = [s;zeros(2^p2-n,1)];     % pading the signal by zeros
n = length(s);
% Padding
sleft = flipud(s(2:n/2+1));
sright = flipud(s(end-n/2:end-1));
x = [sleft; s; sright];
N = length(x);
n1 = length(sleft);
clear xleft xright;
V0 = zeros(n,n);
%%% computes STFT with "Gaussian window" in time domain
ghn0 = @(ci,t) 1/(sqrt(2*pi)*ci)*exp(-t.^2/(2*ci^2));
t = [-n/2+1:n/2-1]/n;  % time sampling index
psih0 = ghn0(ci,t);
x_w = x(n1+ti-n/2+1:n1+ti+n/2-1);
V0 = fft(x_w.*psih0',n); 
F1 = V0(1:n/2)/(n/2);
end

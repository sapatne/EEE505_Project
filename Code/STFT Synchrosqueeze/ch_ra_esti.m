function rk = ch_ra_esti(x,ak,kk,ci,lmd)
% CH_RA_ESTI estimates the chirp rate of a component
%
%  Input:
%   x:  the signal to be analyzed
%   ak: a certain component
%   kk: the time instant
%   ci: the parameter \sigma
%   lmd:   a parameter to determine the length of a window
%
% Output:
%   rk: estimates the chirp rate
%
% Example:
% rk = ch_ra_esti(x,50,128,0.02,1/5);
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
%
x = x(:);
x = x';
arfa = 1/(2*pi)*sqrt(2*log(1/lmd));
n = length(x);
ck = ak;
Ld_h = ceil(2*pi*arfa*ci);           % half duration
tk1 = kk-Ld_h;
tk2 = kk+Ld_h;
s = x;
if tk1<=0
    s = [fliplr(x(1+1:min(abs(tk1)+1+1,n))),x];
    kk = kk+min(abs(tk1),n-1)+1;
end
if tk2>n
    s = [s,fliplr(x(end-(tk2-n+1):end-1))];
end
F1 = [];
for ti = kk-Ld_h:kk+Ld_h
    SF = stft_ti(s,ci,ti);
    F1 = [F1,SF];                                                          % local CWT
end

rt_m = min([ck/(2*Ld_h),(1-ck)/(2*Ld_h),0.25/256]);                           %% consider the maximum chirp rate
rt = -rt_m:1e-4:rt_m;
tt =-Ld_h:Ld_h;   % time index
tp1 = zeros(1,length(rt));
if length(rt)==1
    rk = rt;
else
    for kk = 1:length(rt)
        fkt(kk,:) =rt(kk)*(tt-0)+ck;  %% (Ld_h,ck) is the center of every possible IF
        dkt(kk,:) = mu./fkt(kk,:);       %% the possible ridges cross (Ld_h,ak) on the time-scale plane
        dk_1(kk,:) = round(nv*log2(dkt(kk,:)));
        for dd = 1:2*Ld_h+1;
            tp1(kk) = tp1(kk) + abs(F1(dk_1(kk,dd),dd));
        end
    end
    [ark,rk] = max(tp1);
    rk = rt(rk);  %% chirp rate
end
end


function [up,lo] = boundary_cwt(f,r,arfa,mu,ci)
% calculate the boundaries of multicomponent on the CWT plane;
% input normalized IFs:f, each row- each component
%       normalized chirp rate: r
%       wavelet parameters: mu,ci
% output: upper bound: up
% output: lower bound: lo

% By Dr. Lin Li, 11/15/2016
% Copyright (c) 2017 Xidian University

% This program can not be used for commercialization without the authorization of its author(s). 

% If you use this toolbox for research, please must cite the following papers:
% [1] Lin Li, Haiyan Cai and Qingtang Jiang, "Adaptive synchrosqueezing transform with a time-varying parameter for non-stationary signal separation," preprint, 2018, arXiv:1812.11364.
% [2] Lin Li, Haiyan Cai, Hongxia Han, Qingtang Jiang, and Hongbing Ji, "Adaptive short-time Fourier transform and synchrosqueezing transform for non-stationary signal separation," preprint, 2018, arXiv:1812.11292.
% [3] Haiyan Cai, Qingtang Jiang, Lin Li, and Bruce W. Suter, "Analysis of adaptive short-time Fourier transform-based synchrosqueezing", preprint, 2018, arXiv:1812.11033.

% For any comment or bug report, please send e-mail to lilin@xidian.edu.cn or jiangq@umsl.edu. 





ml = length(r);
n = length(f);
lo = zeros(ml,n);
up = zeros(ml,n);
for ll = 1:ml   
    if abs(r(ll)) <= 1e-5
        lo(ll,:) = (mu-arfa/ci)./f(ll,:);
        up(ll,:) = (mu+arfa/ci)/f(ll,:);
    else
        lo(ll,:) = (-f(ll,:)+sqrt((f(ll,:)).^2-8*pi*arfa*(arfa-mu*ci)*abs(r(ll))))/(4*pi*ci*arfa*abs(r(ll)));
        up(ll,:) = (f(ll,:)-sqrt((f(ll,:)).^2-8*pi*arfa*(arfa+mu*ci)*abs(r(ll))))/(4*pi*ci*arfa*abs(r(ll)));
    end
end
function [Ws,aj] = cwt_mc_ti(s,mu,ci,nv,ti)
% computes the CWT Ws for a fixed time ti
% Input:
%   s: signal
%   gamma: threshold
%   mu,ci: window parameter
% Output:
%   Ws
% By Dr. Lin Li, 11/15/2016
% Copyright (c) 2017 Xidian University
%
%
% This program can not be used for commercialization without the authorization of its author(s). 

% If you use this toolbox for research, please must cite the following paper:
% [1] Lin Li, Haiyan Cai and Qingtang Jiang, "Adaptive synchrosqueezing transform with a time-varying parameter for non-stationary signal separation," preprint, 2018, arXiv:1812.11364.
% [2] Lin Li, Haiyan Cai, Hongxia Han, Qingtang Jiang, and Hongbing Ji, "Adaptive short-time Fourier transform and synchrosqueezing transform for non-stationary signal separation," preprint, 2018, arXiv:1812.11292.
% [3] Haiyan Cai, Qingtang Jiang, Lin Li, and Bruce W. Suter, "Analysis of adaptive short-time Fourier transform-based synchrosqueezing", preprint, 2018, arXiv:1812.11033.

% For any comment or bug report, please send e-mail to lilin@xidian.edu.cn or jiangq@umsl.edu. 


ol = length(s);  % original length
s = s(:);      % Turn into column vector
s = s - mean(s);      %zero mean
n = length(s);     % checking length of signal
p2 = ceil(log2(n));
s = [s;zeros(2^p2-n,1)];     % pading the signal by zeros
n = length(s);
    
L = log2(n);
na = L*nv;    % number of scales
j = [1:na];
m = [0:n-1];   % time index
aj = 2.^(j/nv);    % scales

% Padding, symmetric
sleft = flipud(s(2:n/2+1));
sright = flipud(s(end-n/2:end-1));
x = [sleft; s; sright];
n1 = length(sleft);
clear xleft xright;

N = length(x);
Ws = zeros(na,1);
dWs = Ws;

x = x(:).';
xh = fft(x);

%%% computes CWT with "Morlet"
cs = pi^(-1/4)*ci^(-1/2)*(1+exp(-mu^2*ci^2)-2*exp(-3/4*mu^2*ci^2))^(-1/2);
psihfn = @(w)cs*(2*pi*ci^2)^(1/2)*(exp(-1/2*ci^2*(w-mu).^2)-exp(-1/2*ci^2*(w.^2+mu^2)));

k = 0:(N-1);
xi = zeros(1,N);
xi(1:N/2+1) = 2*pi/N*[0:N/2];
xi(N/2+2:end) = 2*pi/N*[-N/2+1:-1];

t0 = ti+n1;
F1 = 1/N * exp(-1i*2*pi/N*[0:N-1]*t0);
for ai = 1:na
    a = aj(ai);
    psih = psihfn(a*xi);
    dpsih = (i*xi) .* psih;
    xcpsi = F1*((psih .* xh))';
    Ws(ai) = conj(xcpsi); % conj, different from %%%(ifft(psih .* xh))
end

end




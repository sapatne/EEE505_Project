function FSST = regular_pt_adap_sst(s,gamma,sigma)
% REGULAR_PT_ADAP_SST computes the regular phase transform based adaptive SST. 
%
%
% Input
%      s:   the signal to be analyzed
%      gamma:   the threshold for synchrosqueezing
%      sigma:   the time-varying parameter
%
%
% Output
%       FSST:  the regular-PT adaptive SST
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

n = length(s);
FSST = zeros(n/2,n);
for kk = 1:n
    ci = sigma(kk);
    temp = stft_sst(s,gamma,ci);
    FSST(:,kk) = temp(:,kk);
end

end



function SST = stft_sst(s,gamma,ci)
sigma = ci;
n = length(s);
nv = log2(n);
if mod(nv,1)~=0
    warning('The signal is not a power of two, truncation to the next power');
    s = s(1:2^floor(nv));
end
n = length(s);
s = s(:);

bt = 1:n;
ft = 1:n/2;

nb = length(bt);
neta = length(ft);
% Padding
sleft = flipud(s(2:n/2+1));
sright = flipud(s(end-n/2:end-1));
x = [sleft; s; sright];
clear xleft xright;

% Window definition
t = linspace(-0.5,0.5,n);t=t';
g = (1/((2*pi).^(1/2)*sigma))*exp(-1/(2*sigma^2)*t.^2);
gp = -t./sigma^2.* g; % g'
% Initialization
STFT = zeros(neta,nb);
SST = zeros(neta,nb);
omega = zeros(neta,nb);
df = ft(2)-ft(1);
for b=1:nb
	% STFT, window g
    tmp = (fft(x(bt(b):bt(b)+n-1).*g))/n;
    vg = tmp(ft);
    tmp = fft(x(bt(b):bt(b)+n-1).*gp)/n;
    omega(:,b) = (ft-1)'-real(tmp(ft)/2/1i/pi./vg);
    STFT(:,b) = vg.* exp(1i*pi*(ft-1)'); % compensates the tranlation 1/2 of s
end
for b=1:nb
    for eta=1:neta
        if abs(STFT(eta,b))>gamma
            k = 1+round((omega(eta,b)-ft(1)+1)/df);
            if k>=1 && k<=neta
                SST(k,b) = SST(k,b) + STFT(eta,b);% 
            end
        end
    end
end
end

function sig_u = Renyi_entropy(s,sgm,al,t1,mu,nv)
%%% find the time-varying sigma by Renyi entropy
%%% for CWT
% output sig_u:   the time-varying parameter
% input  x:   the signal to be analyzed
%      sgm:   the possible range of the parameter sigma
%       t1:   a parameter to determine the length for integration
%       al:   a constant, al>2
%       mu:   the wavelet paramter
%       nv:   number of voice
%
% Lin Li, October 2017
% Copyright (c) 2017 Xidian University
% 
% 
% This program can not be used for commercialization without the authorization of its author(s). 

% If you use this toolbox for research, please must cite the following paper:
% [1] Lin Li, Haiyan Cai and Qingtang Jiang, "Adaptive synchrosqueezing transform with a time-varying parameter for non-stationary signal separation," preprint, 2018, arXiv:1812.11364.
% [2] Lin Li, Haiyan Cai, Hongxia Han, Qingtang Jiang, and Hongbing Ji, "Adaptive short-time Fourier transform and synchrosqueezing transform for non-stationary signal separation," preprint, 2018, arXiv:1812.11292.
% [3] Haiyan Cai, Qingtang Jiang, Lin Li, and Bruce W. Suter, "Analysis of adaptive short-time Fourier transform-based synchrosqueezing", preprint, 2018, arXiv:1812.11033.

% For any comment or bug report, please send e-mail to lilin@xidian.edu.cn or jiangq@umsl.edu. 

N = length(s);
n_sgm = length(sgm);
Re = zeros(n_sgm,N);
rb2 = zeros(1,N);
rb1 = rb2;

for k = 1:n_sgm
    ci = sgm(k);
    Ws = cwt_mc(s,ci,mu,nv);

    [a,b]=size(Ws);

    alfa=2*al;

    for i=1:t1
        rb1(i)=sum(abs(Ws(:,i)).^alfa);       
        rb2(i)=(sum(abs(Ws(:,i)).^2)).^al;    
    end
    
    for i=N-t1+1:N
        rb1(i)=sum(abs(Ws(:,i)).^alfa);       
        rb2(i)=(sum(abs(Ws(:,i)).^2)).^al;    
    end
    
    for i=t1+1:N-t1
        rb1(i)=sum(sum(abs(Ws(:,(i-t1):(i+t1)))).^alfa);
        rb2(i)=(sum(sum(abs(Ws(:,(i-t1):(i+t1))).^2))).^al;
    end
    for i=1:N
        Re(k,i)=1/(1-al)*log2(rb1(:,i)./rb2(:,i));
    end
end

[x,y] = min(abs(Re),[],1);

sig_u = sgm(y);          %找出每时刻最小entropy值对应的sigma参数值
% sig_u = smooth(sig_u,40,'moving');  %低通滤波平滑
end



function Ws = cwt_mc(s,ci,mu,nv)
% Calculate the CWT    
    Ts = 1;
    Fs = 1;
    s = s(:);          % Turn into column vector
    s = s - mean(s);   % zero mean
    n = length(s);     % checking length of signal
    p2 = ceil(log2(n));
    s = [s;zeros(2^p2-n,1)]; % pading the signal by zeros
    n = length(s);
    
    L = log2(n);
    na = L*nv;         % number of scales
    j = [1:na];
    aj = 2.^(j/nv)*Ts;    % scales
    na = length(aj);
    % Padding, symmetric
    sleft = flipud(s(2:n/2+1));
    sright = flipud(s(end-n/2:end-1));
    x = [sleft; s; sright];
    n1 = length(sleft)-1;
    clear xleft xright;
    N = length(x);
    
    x = x(:).';
    xh = fft(x);
%%% computes CWT with "Morlet"
% cs = pi^(-1/4)*ci^(-1/2)*(1+exp(-mu^2*ci^2)-2*exp(-3/4*mu^2*ci^2))^(-1/2);
% psihfn = @(w)cs*(2*pi*ci^2)^(1/2)*(exp(-1/2*ci^2*(w-mu).^2)-exp(-1/2*ci^2*(w.^2+mu^2)));
    psihfn = @(w) exp(-2*pi^2*ci^2*(w-mu).^2);
    k = 0:(N-1);
    xi = zeros(1,N);
    xi(1:N/2+1) = [0:N/2]/N*Fs;
    xi(N/2+2:end) = [-N/2+1:-1]/N*Fs;
%%%% Compute CWT
    Ws = zeros(na,N);
    for ai = 1:na
        a = aj(ai);
        psih = psihfn(a*xi);
        xcpsi = (ifft(psih .* xh));
        Ws(ai, :) = xcpsi;
    end
    Ws = Ws(:, n1+1:n1+n);
end
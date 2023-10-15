function [Ws Ts1_vt Ts2_vt w2nd aj] = cwt_2nd_order_sst_tv(s,gamma,mu,ci_tv,nv)
% computes the Continuous Wavelet transform: Ws, the 1st-order time-varying synchrosqueezing transform: Ts1_vt.
% the 2nd-order time-varying synchrosqueezing transform: Ts2_vt.   
% and the phase transform w2nd in Eq.(37).
% Input:
%   s: signal
%   gamma: threshold
%   mu: the wavelet parameter
%   ci_tv: the time-varying parameter
%   nv: number of voice
% Output:
%   Ws, Ts1_vt, Ts2_vt
%   w2nd
%
%  By Dr. Lin Li, February 2018
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
s = s - mean(s);   % zero mean
n = length(s);     % checking length of signal
p2 = ceil(log2(n));
s = [s;zeros(2^p2-n,1)];     % pading the signal by zeros
n = length(s);
ci_tv = ci_tv(:);
ci_tv = [ci_tv;zeros(n-ol,1)];

L = log2(n); 
na = L*nv;    % number of scales
j = [1:na];
t = [0:n-1];   % time index
aj = 2.^(j/nv);    % scales

% Padding, symmetric
sleft = flipud(s(2:n/2+1));
sright = flipud(s(end-n/2:end-1));
x = [sleft; s; sright];
n1 = length(sleft);
clear xleft xright;

N = length(x);
W0_t = zeros(na,N);
W11_t = W0_t;
W12_t = W0_t;
W2_t = W0_t;
W0 = zeros(na,ol);
W11 = W0;
W12 = W0;
W2 = W0;

x = x(:).';
xh = fft(x);

%%% compute time-varying CWT with "Morlet's"

psihfn0 = @(w,ci) exp(-2*pi^2*ci^2*(w-mu).^2);  % new definition
psihfn11 = @(w,ci) (2i*pi)*w.*exp(-2*pi^2*ci^2*(w-mu).^2);
psihfn12 = @(w,ci) (-4*pi^2*ci)*(w-mu).^2.*exp(-2*pi^2*ci^2*(w-mu).^2);
psihfn3 = @(w,ci) (-1+4*pi^2*ci^2*(w-mu).^2).*exp(-2*pi^2*ci^2*(w-mu).^2); % \ga 2

xi = zeros(1,N);
xi(1:N/2+1) = 1/N*[0:N/2];
xi(N/2+2:end) = 1/N*[-N/2+1:-1];  % frequency bin

for kk = 1:ol
    ci = ci_tv(kk);
    for ai = 1:na
        a = aj(ai);
        psih0 = psihfn0(a*xi,ci);   
        xcpsi = (ifft(psih0.* xh));
        W0_t(ai, :) = xcpsi;
        
        psih11 = psihfn11(a*xi,ci);   
        xcpsi = (ifft(psih11.* xh));
        W11_t(ai, :) = xcpsi;
        
        psih12 = psihfn12(a*xi,ci);   
        xcpsi = (ifft(psih12.* xh));
        W12_t(ai, :) = xcpsi;    
        
        psih2 = psihfn3(a*xi,ci);   
        xcpsi = (ifft(psih2.* xh));
        W2_t(ai, :) = xcpsi; 
    end
    W0(:,kk) = W0_t(:,n1+kk);
    W11(:,kk) = W11_t(:,n1+kk);
    W12(:,kk) = W12_t(:,n1+kk);
    W2(:,kk) = W2_t(:,n1+kk);
end

%%% compute time-varying SST with "Morlet's" for Eq.(26)
a=repmat(aj,[ol 1]); 
a=a';
cit = ci_tv(:);
cit(ol+1:end) = [];
dcit = ([cit(2:end); 0]-[0;cit(1:end-1)])/2;
dcit(1) = [cit(2)-cit(1)];
dcit(end) = [cit(end)-cit(end-1)];
dcit1 = repmat(dcit,[1,na]);
dcit1 = dcit1';

Ws = W0;
W1 = W11./a + dcit1.*W12;

hWs = (abs(Ws)>gamma).*(1./Ws);
Rs = real(1/(1i*2*pi)*W1.*hWs);% + dcit1./cit .* real(1/(1i*2*pi)*W2.*hWs);    % phase transform
dt = t(2)-t(1);
dT = t(end)-t(1);
fM = 1/(2*dt);
fm = 1/dT;
fs = linspace(fm, fM, na);
dfs = 1/(fs(2)-fs(1));
Ts1 = zeros(length(fs),size(Ws,2));

for b=1:ol
	for ai=1:length(aj)
        if abs(Ws(ai, b)) > gamma
            k = min(max(round(Rs(ai,b)*dfs),1),length(fs));
            Ts1(k, b) = Ts1(k, b) + Ws(ai, b);
        end
    end
end
Ts1_vt = Ts1 * log(2)/nv;

%%% compute time-varying second-order SST (SST2) with "Morlet's" 
%%% by Prof.Jiang (37)
N = length(x);
W0_t = zeros(na,N);
W11_t = W0_t;
W12_t = W0_t;
W2_t = W0_t;
W3_t = W0_t;
W4_t = W0_t;
W5_t = W0_t;
W6_t = W0_t;
W71_t = W0_t;
W72_t = W0_t;
W73_t = W0_t;

W0 = zeros(na,ol);
W11 = W0;
W12 = W0;
W2 = W0;
W3 = W0;
W4 = W0;
W5 = W0;
W6 = W0;
W71 = W0;
W72 = W0;
W73 = W0;

x = x(:).';
xh = fft(x);

psihfn0 = @(w,ci) exp(-2*pi^2*ci^2*(w-mu).^2);  % new definition
psihfn11 = @(w,ci) (2i*pi)*w.*exp(-2*pi^2*ci^2*(w-mu).^2);
psihfn12 = @(w,ci) (-4*pi^2*ci)*(w-mu).^2.*exp(-2*pi^2*ci^2*(w-mu).^2);
psihfn2 = @(w,ci) (2i*pi*ci)*(w-mu).*exp(-2*pi^2*ci^2*(w-mu).^2);
psihfn3 = @(w,ci) (-1+4*pi^2*ci^2*(w-mu).^2).*exp(-2*pi^2*ci^2*(w-mu).^2); 
psihfn4 = @(w,ci) -4*pi^2*ci^2*(w-mu).*w.*exp(-2*pi^2*ci^2*(w-mu).^2);
psihfn5 = @(w,ci) (2i*pi*ci-8i*pi^3*ci^3*(w-mu).^2) .*w.*exp(-2*pi^2*ci^2*(w-mu).^2); 
psihfn6 = @(w,ci) (-16*pi^2*ci^2.*(w-mu).^2 + 12*pi^2*ci^2) .*(w-mu).* w .* exp(-2*pi^2*ci^2*(w-mu).^2);
psihfn71 = @(w,ci) (-8i*pi^3*ci^2).*w.^2.*(w-mu).*exp(-2*pi^2*ci^2*(w-mu).^2);
psihfn72 = @(w,ci) (-8*pi^2*ci).*(w-mu).* w .*exp(-2*pi^2*ci^2*(w-mu).^2);
psihfn73 = @(w,ci) (16*pi^4*ci^3).*(w-mu).^3 .* w .*exp(-2*pi^2*ci^2*(w-mu).^2);

xi = zeros(1,N);
xi(1:N/2+1) = 1/N*[0:N/2];
xi(N/2+2:end) = 1/N*[-N/2+1:-1];

for kk = 1:ol
    ci = ci_tv(kk);
    for ai = 1:na
        a = aj(ai);
        
        psih0 = psihfn0(a*xi,ci);   
        xcpsi = (ifft(psih0.* xh));
        W0_t(ai, :) = xcpsi;
        
        psih11 = psihfn11(a*xi,ci);   
        xcpsi = (ifft(psih11.* xh));
        W11_t(ai, :) = xcpsi;
        
        psih12 = psihfn12(a*xi,ci);   
        xcpsi = (ifft(psih12.* xh));
        W12_t(ai, :) = xcpsi;    
        
        psih2 = psihfn2(a*xi,ci);   
        xcpsi = (ifft(psih2.* xh));
        W2_t(ai, :) = xcpsi; 

        psih3 = psihfn3(a*xi,ci);   
        xcpsi = (ifft(psih3.* xh));
        W3_t(ai, :) = xcpsi; 
        
        psih4 = psihfn4(a*xi,ci);   
        xcpsi = (ifft(psih4.* xh));
        W4_t(ai, :) = xcpsi;         
        
        psih5 = psihfn5(a*xi,ci);   
        xcpsi = (ifft(psih5.* xh));
        W5_t(ai, :) = xcpsi; 
        
        psih6 = psihfn6(a*xi,ci);   
        xcpsi = (ifft(psih6.* xh));
        W6_t(ai, :) = xcpsi;
        
        psih71 = psihfn71(a*xi,ci);   
        xcpsi = (ifft(psih71.* xh));
        W71_t(ai, :) = xcpsi;
        
        psih72 = psihfn72(a*xi,ci);   
        xcpsi = (ifft(psih72.* xh));
        W72_t(ai, :) = xcpsi;
        
        psih73 = psihfn73(a*xi,ci);   
        xcpsi = (ifft(psih73.* xh));
        W73_t(ai, :) = xcpsi;
        
    end
    W0(:,kk) = W0_t(:,n1+kk);
    W11(:,kk) = W11_t(:,n1+kk);
    W12(:,kk) = W12_t(:,n1+kk);
    W2(:,kk) = W2_t(:,n1+kk);
    W3(:,kk) = W3_t(:,n1+kk);
    W4(:,kk) = W4_t(:,n1+kk);
    W5(:,kk) = W5_t(:,n1+kk);
    W6(:,kk) = W6_t(:,n1+kk);
    W71(:,kk) = W71_t(:,n1+kk);
    W72(:,kk) = W72_t(:,n1+kk);
    W73(:,kk) = W73_t(:,n1+kk);
end

% phase transform
a=repmat(aj,[ol 1]); 
a=a';

W1 = W11./a + dcit1.*W12;
W4 = W4./a;
W5 =  W5./a;
W6 = W6./a;
W7 = 1./(a.^2) .* W71 + dcit1./a .*(W72+W73);

tp1 = W2./W0 + a.*(W5.*W0-W2.*W4)./W0.^2;
tp2 = (W7.*W0-W1.*W4)./W0.^2;
tp3 = (W6.*W0-W3.*W4)./W0.^2;

w2nd = zeros(na,ol); 
R0 = zeros(na,ol);
for kk = 2:ol-1
    for ai = 2:na-1
        if abs(Ws(ai,kk)) >= gamma
            R0(ai,kk) = 1./tp1(ai,kk) * (tp2(ai,kk)+dcit(kk)/cit(kk)*tp3(ai,kk));
            if abs(R0(ai,kk)) > 0
                w2nd(ai,kk) = real(W1(ai,kk)/(2i*pi*W0(ai,kk))) - a(ai,kk)*real(W2(ai,kk)/(2i*pi*W0(ai,kk)) .* R0(ai,kk)) + dcit(kk)/cit(kk) * real(W3(ai,kk)/(2i*pi*W0(ai,kk)));
            else
                w2nd(ai,kk) = real(W1(ai,kk)/(2i*pi*W0(ai,kk))) + dcit(kk)/cit(kk) * real(W3(ai,kk)/(2i*pi*W0(ai,kk)));
            end
        end
    end
end

dT = t(end)-t(1);
fM = 1/(2*dt);
fm = 1/dT;
fs = linspace(fm, fM, na);
dfs = 1/(fs(2)-fs(1));

Ts2 = zeros(length(fs),size(Ws,2));
for b=1:ol
    for ai=1:length(aj)
        if abs(W0(ai, b)) > gamma
            k = min(max(round(w2nd(ai,b)*dfs),1),length(fs));
            Ts2(k, b) = Ts2(k, b) + W0(ai, b);
        end
    end
end
Ts2_vt = Ts2 * log(2)/nv;
end




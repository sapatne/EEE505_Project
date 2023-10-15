function [Vs_tv Ts1_tv Ts2_tv] = adap_stft_sst(s,gamma,sigma_tv)
% ADAP_STFT_SST computes the adaptive STFT Vs_tv, the adaptive FSST
% Ts_tv and adaptive 2nd-order FSST with time-varying parameter. 
%
% Input:
%   s: signal
%   gamma: threshold
%   sigma_vt: time-varying parameter for Gaussian window 
% 
% Output:
%   Vs_tv: adaptive STFT
%   Ts_tv: adaptive SST   
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
% N = n; % no padding;
% n1 = 0;% no padding;
% x = s; % no padding;

cit = sigma_tv(:);
cit1 = repmat(cit,[1,n]);
cit1 = cit1';
dcit = ([cit(2:end); 0]-[0;cit(1:end-1)])/2;
dcit(1) = [cit(2)-cit(1)];
dcit(end) = [cit(end)-cit(end-1)];
dcit1 = repmat(dcit,[1,n]);
dcit1 = dcit1';       % the derivative of \sigma(t)     

%%%% compute different types of time-varying STFT
V0 = zeros(N,n);  % g
V1 = zeros(N,n);  % g1
V2 = zeros(N,n);  % g2
V3 = zeros(N,n);  % derivative of V0 with respect to t
V4 = zeros(N,n);  % derivative of V0 with respect to \xi 
V5 = zeros(N,n);  % derivative of V1 with respect to \xi
V6 = zeros(N,n);  % derivative of V2 with respect to \xi
V7 = zeros(N,n);  % second order derivative of V0 with respect to t and \xi

V0_t = zeros(N,N);
V1_t = zeros(N,N);
V2_t = zeros(N,N);
V31_t = zeros(N,N);
V32_t = zeros(N,N);
V4_t = zeros(N,N);
V5_t = zeros(N,N);
V6_t = zeros(N,N);
V71_t = zeros(N,N);
V72_t = zeros(N,N);
V73_t = zeros(N,N);

%%% FFT of x(t)
x = x(:).';
xh = fft(x);
% sh = fftshift(sh);
xi = zeros(1,N);
xi(1:N/2+1) = 1/N*[0:N/2];
xi(N/2+2:end) = 1/N*[-N/2+1:-1];


%%% computes STFT with "Gaussian window" in the frequency domain
ghn0 = @(ci,xi) exp(-2*pi^2*ci^2*xi.^2);
ghn1 = @(ci,xi) -2i*pi*ci*xi.* exp(-2*pi^2*ci^2*xi.^2);
ghn2 = @(ci,xi) (4*pi^2*ci^2*xi.^2-1).* exp(-2*pi^2*ci^2*xi.^2);
ghn3 = @(ci,xi) -4*pi^2*ci*(xi).^2 .* exp(-2*pi^2*ci^2*xi.^2);
ghn4 = @(ci,xi) -4*pi^2*ci^2*(xi).* exp(-2*pi^2*ci^2*xi.^2);
ghn5 = @(ci,xi) 8i*pi^3*ci^3*(xi).^2 .* exp(-2*pi^2*ci^2*xi.^2);
ghn6 = @(ci,xi) -4*pi^2*ci^2*(xi).* (4*pi^2*ci^2*xi.^2-1).* exp(-2*pi^2*ci^2*xi.^2);
ghn71 = @(ci,xi) -8i*pi^3*ci^2*(xi).* exp(-2*pi^2*ci^2*xi.^2);
ghn72 = @(ci,xi) -8*pi^2*ci*(xi).* exp(-2*pi^2*ci^2*xi.^2);
ghn73 = @(ci,xi) 16*pi^4*ci^3*(xi).^3 .* exp(-2*pi^2*ci^2*xi.^2);

%%%
for bb = 1:ol
    ci = sigma_tv(bb);
    for kk = 1:N
        eta = xi(kk);  
        
        psih0 = ghn0(ci,eta-xi);
        xcpsi = (ifft((psih0.* xh)));
        V0_t(kk,:) = xcpsi;
        
        psih1 = ghn1(ci,eta-xi);
        xcpsi = (ifft((psih1.* xh)));
        V1_t(kk,:) = xcpsi;
        
        psih2 = ghn2(ci,eta-xi);
        xcpsi = (ifft((psih2.* xh)));
        V2_t(kk,:) = xcpsi;
        
        psih31 = 2i*pi*xi.*ghn0(ci,eta-xi);
        xcpsi = (ifft((psih31.* xh)));
        V31_t(kk,:) = xcpsi;
        
        psih32 = ghn3(ci,eta-xi);
        xcpsi = dcit(bb)*(ifft((psih32.* xh)));
        V32_t(kk,:) = xcpsi;
        
        psih4 = ghn4(ci,eta-xi);
        xcpsi = (ifft((psih4.* xh)));
        V4_t(kk,:) = xcpsi;
        
        psih5 = ghn5(ci,eta-xi);
        xcpsi = (ifft((psih5.* xh)));
        V5_t(kk,:) = xcpsi-2i*pi*ci*V0_t(kk,:);
        
        psih6 = ghn6(ci,eta-xi);
        xcpsi = (ifft((psih6.* xh)));
        V6_t(kk,:) = xcpsi+4i*pi*ci*V1_t(kk,:);
        
        psih71 = xi.* ghn71(ci,eta-xi);
        xcpsi = (ifft((psih71.* xh)));
        V71_t(kk,:) = xcpsi;

        psih72 = ghn72(ci,eta-xi);
        xcpsi = dcit(bb) * (ifft((psih72.* xh)));
        V72_t(kk,:) = xcpsi;

        psih73 = ghn73(ci,eta-xi);
        xcpsi = dcit(bb) * (ifft((psih73.* xh)));
        V73_t(kk,:) = xcpsi;     
    end
    V0(:,bb) = V0_t(:,bb+n1);
    V1(:,bb) = V1_t(:,bb+n1);
    V2(:,bb) = V2_t(:,bb+n1);   
    V3(:,bb) = V31_t(:,bb+n1) + V32_t(:,bb+n1);
    V4(:,bb) = V4_t(:,bb+n1);   
    V5(:,bb) = V5_t(:,bb+n1);   
    V6(:,bb) = V6_t(:,bb+n1);
    V7(:,bb) = V71_t(:,bb+n1) + V72_t(:,bb+n1) + V73_t(:,bb+n1);
end

% for real signals:
V0(n+1:end,:) = [];
V1(n+1:end,:) = [];
V2(n+1:end,:) = [];
V3(n+1:end,:) = [];
V4(n+1:end,:) = [];
V5(n+1:end,:) = [];
V6(n+1:end,:) = [];
V7(n+1:end,:) = [];

% time-varying STFT
Vs_tv = V0/(n/2);

% adative FSST (first order)
Ts1_tv = zeros(n/2,ol);
w1st = real(V3./(2i*pi*V0)) + dcit1./cit1 .* real(V2./(2i*pi*V0)); % Eq.(20)
w1st = real(w1st) .* (abs(V0)>gamma);  % real Jan 1, 2019
fs = [0:n/2-1]/n;
dfs = n;

for bb=1:ol
	for fi=1:n
        if abs(V0(fi, bb)) > gamma
            kk = min(max(round(w1st(fi,bb)*dfs),1),length(fs));
            Ts1_tv(kk, bb) = Ts1_tv(kk, bb) + V0(fi, bb) * cit(bb);
        end
    end
end
Ts1_tv = Ts1_tv/(dfs/2);


%%% second-order adaptive SST
Ts2_tv = zeros(n/2,ol);
w2nd = zeros(n/2,ol); % Eq.(30) 
tp1 = (V5.*V0-V1.*V4)./V0.^2;
tp2 = (V7.*V0-V3.*V4)./V0.^2;
tp3 = (V6.*V0-V2.*V4)./V0.^2;
% P0 = 1./tp1 .* (tp2 + dcit1./cit1 * tp3);    % Eq.(31), forgot a '.' here!!!!!!
R0 = real(V3./(2i*pi*V0));% .* (abs(V0)>gamma);
R1 = dcit1./cit1 .* real(V2./(2i*pi*V0));%.* (abs(V0)>gamma);
R2 = V1./(2i*pi*V0);%.* (abs(V0)>gamma);

for bb = 1:ol
    for fi = 1:n
        if tp1(fi,bb) >= gamma
            P0 = (tp2(fi,bb) + dcit(bb)/cit(bb) * tp3(fi,bb))/tp1(fi,bb);    % Eq.(31)
            w2nd(fi,bb) = R0(fi,bb) + R1(fi,bb) - real(R2(fi,bb) * P0);
        else
            w2nd(fi,bb) = R0(fi,bb) + R1(fi,bb);
        end
    end
end

w2nd = real(w2nd); % real, Jan 1, 2019
fs = [0:n/2-1]/n;
dfs = n;
for bb=1:ol
	for fi=1:n
        if abs(V0(fi, bb)) > gamma
            kk = min(max(round(w2nd(fi,bb)*dfs),1),length(fs));
            Ts2_tv(kk, bb) = Ts2_tv(kk, bb) + V0(fi, bb) * cit(bb);
        end
    end
end
Ts2_tv = Ts2_tv/(dfs/2);




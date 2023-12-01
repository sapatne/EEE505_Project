% ASST toolbox v1.0
% this is an example for a two-component signal
% Lin Li, October 2017
% Copyright (c) 2017 Xidian University
% This program can not be used for commercialization without the authorization of its author(s). 

% If you use this toolbox for research, please must cite the following papers:
% [1] Lin Li, Haiyan Cai and Qingtang Jiang, "Adaptive synchrosqueezing transform with a time-varying parameter for non-stationary signal separation," preprint, 2018, arXiv:1812.11364.
% [2] Lin Li, Haiyan Cai, Hongxia Han, Qingtang Jiang, and Hongbing Ji, "Adaptive short-time Fourier transform and synchrosqueezing transform for non-stationary signal separation," preprint, 2018, arXiv:1812.11292.
% [3] Haiyan Cai, Qingtang Jiang, Lin Li, and Bruce W. Suter, "Analysis of adaptive short-time Fourier transform-based synchrosqueezing", preprint, 2018, arXiv:1812.11033.

% For any comment or bug report, please send e-mail to lilin@xidian.edu.cn or jiangq@umsl.edu. 


clear;
close all;
clc;
%%%%%%% 2018.1.14
N = 256;

Fs = 256;
Ts = 1/Fs;
t = [0:N-1]/Fs; % s
fi = t.*Fs;    % Hz  

% r1 = 50;r2 =64; c1 = 12;c2 =34;
% 
% x1 = cos(2*pi*(c1*t+r1/2*t.^2));
% x2 = cos(2*pi*(c2*t+r2/2*t.^2));
% f1 = c1+r1*t;
% f2 = c2+r2*t;
% s = x1+x2;

[s, f, r, c] = signal_gen('n');
r1 = r(1);
r2 = r(2);
c1 = c(1);
c2 = c(2);
f1 = f(1, :);
f2 = f(2, :);


figure;
plot(t,s,'k-','linewidth',1);
title("Signal in Time domain")

figure;
plot(t,f1,'r-','linewidth',2);
hold on;
plot(t,f2,'b-','linewidth',2);
title("Instantenous frequencies")

% figure;
% plot(fi,abs(fft(hilbert(s))),'k-','linewidth',1);
% axis([0 0.5 0 75]);
% xlabel('Frquency (Hz)');
% ylabel('Amplitude');

gamma = 0.01;
mu = 1;
%ci = 2;
ci = 0.5; 
nv = 32;
ol = length(s);    % original length
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
% aj = 1:0.1:30;
na = length(aj);

% Padding, symmetric
sleft = flipud(s(2:n/2+1));
sright = flipud(s(end-n/2:end-1));
x = [sleft; s; sright];
n1 = length(sleft);
clear xleft xright;

N = length(x);

x = x(:).';
xh = fft(x);
%% computes CWT with "Morlet"
psihfn = @(w) exp(-2*pi^2*ci^2*(w-mu).^2);
k = 0:(N-1);
xi = zeros(1,N);
xi(1:N/2+1) = [0:N/2]/N*Fs;
xi(N/2+2:end) = [-N/2+1:-1]/N*Fs;

%% Compute CWT
Ws = zeros(na,N);
dWs = Ws;
for ai = 1:na
    a = aj(ai);
    psih = psihfn(a*xi);
    dpsih = (i*xi) .* psih;

    xcpsi = (ifft(psih .* xh));
    Ws(ai, :) = xcpsi;
    dxcpsi = (ifft(dpsih .* xh));
    dWs(ai, :) = dxcpsi;
end
Ws = Ws(:, n1+1:n1+n);
dWs = dWs(:, n1+1:n1+n);

lmd = 1/5;    %%duration
%lmd = 1/10; 
arfa = 1/(2*pi)*sqrt(2*log(1/lmd)); 


figure;
imagesc(t,aj,abs(Ws));
title("CWT with Morlet")

%% Calculate sigma_1 and sigma_2

% the lower edge for x1
l1 = (-f1+sqrt(f1.^2-8*pi*arfa*(arfa-mu*ci)*abs(r1)))/(4*pi*ci*arfa*abs(r1));
% the upper edge for x2
u2 = (f2-sqrt(f2.^2-8*pi*arfa*(arfa+mu*ci)*abs(r2)))/(4*pi*ci*arfa*abs(r2));
jl1 = nv*log2(l1/Ts);
ju2 = nv*log2(u2/Ts);
hold on;
plot(t,jl1/Fs,'k-','Linewidth',2);
plot(t,ju2/Fs,'b--','Linewidth',2);

% Sigma_1
ci_1 = arfa/mu * (f1+f2)./(abs(f2-f1));

% Sigma_2
r1=abs(r1);r2=abs(r2); 
ar1 = 2*pi*arfa*mu*(r1+r2).^2;
bt1 = (r1*f2+r2*f1).*(f2-f1)+4*pi*arfa^2*(r2^2-r1^2);
gm1 = arfa/mu*((r1*f2+r2*f1).*(f2+f1)+2*pi*arfa^2*(r2-r1)^2); 
dt1=(r1*f2+r2*f1).^2.*((f2-f1).^2-16*pi*arfa^2*(r1+r2));
ci_2 = max(arfa/mu,real((bt1-sqrt(dt1))/(2*ar1)));

ci_2up = (bt1+sqrt(dt1))/(2*ar1);
%%

figure;
plot(t,ci_1,'k-','linewidth',2);
hold on;
plot(t,ci_2,'r-','linewidth',2);

%% Compute SST with sigma = 1
[Ws, Rs, STs, STs2, aj] = cwt_sst_mc(s,gamma,mu*2*pi,ci,nv);
figure;
subplot(431);
imagesc(t,aj*Ts,abs(Ws));
title("CWT 2")

subplot(432);
imagesc(t,[0:n/2-1],abs(STs));
title("SST of CWT")
subplot(433);
imagesc(t,[0:n/2-1],abs(STs2));
title("SST 2 of CWT")

%% Compute Time-varying parameter SST with sigma_1
for kk = 1:n
    ci = ci_1(kk);
    [Ws Rs STs aj] = cwt_sst_mc_ti(s,gamma,mu*2*pi,ci,nv,kk-1);
    Ws_1(:,kk) = Ws;
    STs_1(:,kk) = STs;
end
subplot(434);
imagesc(t,aj*Ts,abs(Ws_1));
title("CWT 3")

l1_p1 = (mu-arfa./ci_1)./f1;           % signals are approximated by harmanics
u2_p1 = (mu+arfa./ci_1)./f2;
jl1_p1 = nv*log2(l1_p1/Ts);
ju2_p1 = nv*log2(u2_p1/Ts);

hold on;
plot(t,jl1_p1/Fs,'k-','Linewidth',2);
plot(t,ju2_p1/Fs,'b--','Linewidth',2);

subplot(435);
imagesc(abs(STs_1));
title("Time Varying SST of CWT")

% Time-varying SST with sigma_2
for kk = 1:n
    ci = ci_2(kk);
    [Ws Rs STs aj] = cwt_sst_mc_ti(s,gamma,mu*2*pi,ci,nv,kk-1);
    Ws_2(:,kk) = Ws;
    STs_2(:,kk) = STs;
end
% subplot(436);
% imagesc(t,aj*Ts,abs(Ws_2));
% title("CWT 4")

l1_p = (-f1+sqrt(f1.^2-8*pi*arfa*(arfa-mu*ci_2)*abs(r1)))./(4*pi*ci_2*arfa*abs(r1));
% the upper edge for x2
u2_p = (f2-sqrt(f2.^2-8*pi*arfa*(arfa+mu*ci_2)*abs(r2)))./(4*pi*ci_2*arfa*abs(r2));
jl1_p = nv*log2(l1_p/Ts);
ju2_p = nv*log2(u2_p/Ts);

hold on;
plot(t,jl1_p/Fs,'k-','Linewidth',2);
plot(t,ju2_p/Fs,'b--','Linewidth',2);
legend('\it{l}\rm{_1}','\it{u}\rm{_2}');
set(legend,'FontName','Times New Roman')


subplot(436);
imagesc(abs(STs_2));
title("Time Vary SST of CWT 2")


%% compute the CWT, conventional SST, conventional 2nd-order SST.
ci_tv = ones(1, N);
[Ws, Ts1_vt, Ts2_vt, w2nd, aj] = cwt_2nd_order_sst_tv(s,gamma,mu,ci_tv,nv);

subplot(437);
imagesc(t, aj*Ts, abs(Ws))
title("CWT for \sigma=1")
axis xy;

subplot(438);
imagesc(t, [0:n/2-1], abs(Ts1_vt));
title("1st order SST for \sigma=1")
axis xy;

subplot(439);
imagesc(t, [0:n/2-1], abs(Ts2_vt));
title("2nd order SST for \sigma=1")
axis xy;

%% compute the adaptive CWT, adaptive SST, adaptive 2nd-order SST.
ci_tv = ci_2;
[Ws, Ts1_vt, Ts2_vt, w2nd, aj] = cwt_2nd_order_sst_tv(s,gamma,mu,ci_tv,nv);

% subplot(4310);
imagesc(t, aj*Ts, abs(Ws))
title("Adaptive CWT with \sigma_2")
axis xy;

% subplot(4311);
imagesc(t, [0:n/2-1], abs(Ts1_vt));
title("Adaptive 1st order SST for \sigma_2")
axis xy;

% subplot(4312);
imagesc(t, [0:n/2-1], abs(Ts2_vt));
title("Adaptive 2nd order SST for \sigma_2")
axis xy;
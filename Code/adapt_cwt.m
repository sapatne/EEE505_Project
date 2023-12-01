
clear;
close all;
clc;

Fs = 256;
Ts = 1/Fs;
N = Fs;
t = linspace(0, 1-1/N, N); % s
fi = t.*Fs;    % Hz  

%%
% c1 = 25;
% r1 = 100;
% 
% c2 = 64;
% r2 = 12; 
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

% figure;
% plot(t,s,'k-','linewidth',1);
% title("Signal in Time domain")

figure;
plot(t,f1,'r-','linewidth',2);
hold on;
plot(t,f2,'b-','linewidth',2);
title("Instantenous frequencies")

%% Initialize params

gamma = 0.01;
mu = 1;
ci = 0.1; 
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

%% Calculate sigma1 & sigma2

lmd = 1/5;         % duration
arfa = 1/(2*pi)*sqrt(2*log(1/lmd)); 

% Sigma_1
ci_1 = arfa/mu * (f1+f2)./(f2-f1);

% Sigma_2
% r1=abs(r1);r2=abs(r2); 
ar1 = 2*pi*arfa*mu*(r1+r2).^2;
bt1 = (r1*f2+r2*f1).*(f2-f1)+4*pi*arfa^2*(r2^2-r1^2);
gm1 = arfa/mu*((r1*f2+r2*f1).*(f2+f1)+2*pi*arfa^2*(r2-r1)^2); 
dt1=(r1*f2+r2*f1).^2.*((f2-f1).^2-16*pi*arfa^2*(r1+r2));
ci_2 = max(arfa/mu,real((bt1-sqrt(dt1))/(2*ar1)));

ci_2up = (bt1+sqrt(dt1))/(2*ar1);

%% compute the CWT, conventional SST, conventional 2nd-order SST.

ci_tv = ones(1, N);
[Ws, Ts1_vt, Ts2_vt, w2nd, aj] = cwt_2nd_order_sst_tv(s,gamma,mu,ci_tv,nv);

subplot(231)
imagesc(t, aj*Ts, abs(Ws))
title("CWT for fixed \sigma=1")
axis xy;

% figure;
subplot(232)
imagesc(t, [0:n/2-1], abs(Ts1_vt));
title("1st order WSST for fixed \sigma=1")
axis xy;

% figure;
subplot(233)
imagesc(t, [0:n/2-1], abs(Ts2_vt));
title("2nd order WSST for fixed \sigma=1")
axis xy;

%% compute the adaptive CWT, adaptive SST, adaptive 2nd-order SST.

ci_tv = ci_2;
[Ws, Ts1_vt, Ts2_vt, w2nd, aj] = cwt_2nd_order_sst_tv(s,gamma,mu,ci_tv,nv);

% figure;
subplot(234)
imagesc(t, aj*Ts, abs(Ws))
title("CWT with adaptive \sigma")
axis xy;

subplot(235)
imagesc(t, [0:n/2-1], abs(Ts1_vt));
title("1st order WSST for adaptive \sigma")
axis xy;

subplot(236)
imagesc(t, [0:n/2-1], abs(Ts2_vt));
title("2nd order WSST for adaptive \sigma")
axis xy;

% adaptive FSST toolbox v1.0
% this is an example for the two-component LFM signal

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




clear;
close all;
clc;

load Data\LinearChirp.mat signal
Fs = 2048;
N = length(signal);
t = linspace(0,1,N);

c1 = 1000;
c2 = 0;
r1 = (125-1000);
r2 = 0; % the one for adaptive wavelet-based SST

f1 = c1+r1*t;
f2 = 0;
s = signal;

gamma = 0.001;

lmd = 1/5.2;    %%duration
arfa = 1/(2*pi)*sqrt(2*log(1/lmd)); 
ci_1 = arfa./((f2-f1)/N);        % Eq.(47)
tic
Ak = 2*pi*arfa*(abs(r1/(N^2)) + abs (r2/(N^2)));  % normalized
if length(Ak)==1;
    Ak = ones(1,N)*Ak;
end
Bk = (f2-f1)/N;
Ck = sqrt(Bk.^2-8*arfa*Ak);
for bb = 1:N
    if Ak(bb)~=0 && real(Ck(bb)) == Ck(bb)
        ci_2(bb) = max(1,(Bk(bb)-Ck(bb))/(2*Ak(bb)));
    else
        ci_2(bb) = ci_1(bb);
    end
end
   
ci_1 = ci_1/N;
ci_2 = ci_2/N;

%%%% the Renyi Entropy method
sgm_1 = 0.005;
d_sgm = 0.001;
sgm = sgm_1:d_sgm:0.1;
al = 2.5;
t1 = 4;
[sgm_u,sgm_R1,sgm_R2] = Renyi_STFT_SST(s,sgm,al,t1);% Entropy
ci_tv2=smooth(sgm_R1,20,'rlowess');                          % \sigma_{Re}
ci_tv3=smooth(sgm_R2,20,'rlowess');                          % \sigma_{Re2}

gamma1 = 0.3;
ci_est = tv_para_fast(s,sgm_1,sgm_u,d_sgm,lmd,gamma1);  
ci_tv1=smooth(ci_est,20,'rlowess'); 


figure;
plot(t,ci_1,'b--','linewidth',2);
hold on;
plot(t,ci_2,'k--','linewidth',2);
plot(t,sgm_u,'m-','linewidth',2);
plot(t,ci_tv1,'r-','linewidth',2);
plot(t,ci_tv2,'g-','linewidth',2);
plot(t,ci_tv3,'c-','linewidth',2);
legend('\sigma_1(t)', '\sigma_2(t)','\sigma_u(t)','\sigma_{est}(t)','\sigma_{Re}(t)','\sigma_{Re2}(t)','Location','north')
axis([0 1 0 0.11]);
set(gca,'FontSize',20)
xlabel('Time (s)','FontSize',20);
ylabel('\sigma','FontSize',20);
title('Various time-varying parameters of \sigma(t)');

sigma_tv = ci_2;
[Vs_tv Ts1_tv Ts2_tv] = adap_stft_sst(s,gamma,sigma_tv*Fs);
figure;
imageSQ(t,[0:N/2-1],abs(Vs_tv));
axis xy;
xlabel('Time (s)','FontSize',20);
ylabel('Frequency (Hz)','FontSize',20);
set(gca,'FontSize',20);
title('Adaptive STFT with \sigma_2(t)');
figure;
imageSQ(t,[0:N/2-1],abs(Ts1_tv));
axis xy;
xlabel('Time (s)','FontSize',20);
ylabel('Frequency (Hz)','FontSize',20);
set(gca,'FontSize',20);
title('Adaptive FSST with \sigma_2(t)');
figure;
imageSQ(t,[0:N/2-1],abs(Ts2_tv));
axis xy;
xlabel('Time (s)','FontSize',20);
ylabel('Frequency (Hz)','FontSize',20);
set(gca,'FontSize',20);
title('Adaptive 2nd-order FSST with \sigma_2(t)');

sigma_tv = ci_tv1;
[Vs_tv Ts1_tv Ts2_tv] = adap_stft_sst(s,gamma,sigma_tv*Fs);
figure;
imageSQ(t,[0:N/2-1],abs(Vs_tv));
% axis([0 1 0 40]);
axis xy;
xlabel('Time (s)','FontSize',20);
ylabel('Frequency (Hz)','FontSize',20);
set(gca,'FontSize',20);
title('Adaptive STFT with the estimated \sigma(t)');
figure;
imageSQ(t,[0:N/2-1],abs(Ts1_tv));
axis xy;
xlabel('Time (s)','FontSize',20);
ylabel('Frequency (Hz)','FontSize',20);
set(gca,'FontSize',20);
title('Adaptive FSST with the estimated \sigma(t)');
figure;
imageSQ(t,[0:N/2-1],abs(Ts2_tv));
axis xy;
xlabel('Time (s)','FontSize',20);
ylabel('Frequency (Hz)','FontSize',20);
set(gca,'FontSize',20);
title('Adaptive 2nd-order FSST with the estimated \sigma(t)');



% regular-PT adaptive SST with Re1
FSST_rpt = regular_pt_adap_sst(s,gamma,ci_tv2); 
figure;
imageSQ(t,[0:N/2-1],abs(FSST_rpt));
axis xy;
xlabel('Time (s)','FontSize',20);
ylabel('Frequency (Hz)','FontSize',20);
set(gca,'FontSize',20);
title('Regular-PT adaptive SST with \sigma_{Re1}(t)');


% regular-PT adaptive 2nd-order SST with Re2
FSST2_rpt = regular_pt_adap_sst2(s,gamma,ci_tv3);
figure;
imageSQ(t,[0:N/2-1],abs(FSST2_rpt));
axis xy;
xlabel('Time (s)','FontSize',20);
ylabel('Frequency (Hz)','FontSize',20);
set(gca,'FontSize',20);
title('2nd-order regular-PT adaptive SST with \sigma_{Re2}(t)');


[STFT SST VSST] = stft_sst2(s,gamma,0.02);
figure;
imageSQ(t,[0:N/2-1],abs(SST));
axis xy;
xlabel('Time (s)','FontSize',20);
ylabel('Frequency (Hz)','FontSize',20);
set(gca,'FontSize',20);
title('Conventional SST when \sigma=0.02');
figure;
imageSQ(t,[0:N/2-1],abs(VSST));
axis xy;
xlabel('Time (s)','FontSize',20);
ylabel('Frequency (Hz)','FontSize',20);
set(gca,'FontSize',20);
title('Conventional 2nd-order SST when \sigma=0.02');

toc

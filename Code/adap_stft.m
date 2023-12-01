% clear;
% close all;
% clc;

N = 256;
Ts = 1/N;

t = [0:N-1]/N; % s
fi = t.*N;    % Hz  

% c1 = 12;s
% c2 =64;
% 
% r1 = 50;
% r2 =-34; 
% 
% x1 = cos(2*pi*(c1*t+r1/2*t.^2));
% x2 = cos(2*pi*(c2*t+r2/2*t.^2));
% f1 = c1+r1*t;
% f2 = c2+r2*t;
% s = x1+x2;

[s, f, r, c] = signal_gen();
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

%% initial sigma1 & sigma2
gamma = 0.001;

lmd = 1/5.2;    %%duration
arfa = 1/(2*pi)*sqrt(2*log(1/lmd)); 
ci_1 = arfa./((f2-f1)/N);        % Eq.(47)
tic
Ak = 2*pi*arfa*(abs(r1/(N^2)) + abs (r2/(N^2)));  % normalized
if length(Ak)==1
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

%% the Renyi Entropy method
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

%% conventional stft and fsst
[STFT, SST, VSST] = stft_sst2(s,gamma,0.02);

subplot(231);
imagesc(t, fi, abs(STFT))
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Conventional STFT when \sigma=0.02');


subplot(232);
imagesc(t,fi,abs(SST));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Conventional SST when \sigma=0.02');

subplot(233);
imagesc(t,fi,abs(VSST));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Conventional 2nd-order SST when \sigma=0.02');


%% estimated sigma for fsst
sigma_tv = ci_tv1;
[Vs_tv Ts1_tv Ts2_tv] = adap_stft_sst(s,gamma,sigma_tv*N);
subplot(234);
imagesc(t,fi,abs(Vs_tv));
% axis([0 1 0 40]);
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Adaptive STFT with the estimated \sigma(t)');

subplot(235);
imagesc(t,fi,abs(Ts1_tv));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Adaptive FSST with the estimated \sigma(t)');

subplot(236);
imagesc(t,fi,abs(Ts2_tv));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Adaptive 2nd-order FSST with the estimated \sigma(t)');


function [sig_u,sig_sst1,sig_sst2] = Renyi_STFT_SST(s,sgm,al,t1)
% Renyi_STFT_SST  finds the time-varying sigma by Renyi entropy
%
%
% input  s:   the signal to be analyzed
%      sgm:   the possible range of the parameter sigma
%       t1:   a parameter to determine the length for integration
%       al:   a constant, al>2
%
%
% output sig_u:  the time-varying parameter based on the  Renyi entropy of STFT
%        sig_sst1:   the time-varying parameter based on the  Renyi entropy of SST
%        sig_sst2:   the time-varying parameter based on the  Renyi entropy of SST2
%
% Example:
% sgm_1 = 0.005, d_sgm = 0.001;
% sgm = sgm_1:d_sgm:0.1;
% al = 2.5, t1 = 4;
% [sgm_u,sgm_R1,sgm_R2] = Renyi_entropy_STFT_SST(s,sgm,al,t1);% Entropy
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



%check length
if (nargin~=4)
    error('Four input parameters are required')
end

N = length(s);
n_sgm = length(sgm);
Re0 = zeros(n_sgm,N);
Re1 = zeros(n_sgm,N);
Re2 = zeros(n_sgm,N);
rb2 = zeros(1,N);
rb1 = rb2;
gamma = 0.001;
alfa=2*al;

for k = 1:n_sgm
    ci = sgm(k);
    
    [STFT SST SST2] = stft_sst2(s,gamma,ci);

    Ds = abs(STFT);
    Ds = [fliplr(Ds(:,1:t1)) Ds fliplr(Ds(:,end-t1+1:end))]; %padding
    for i=t1+1:N-t1
        rb1(i-t1)=sum(sum(Ds(:,(i-t1):(i+t1)).^alfa));
        rb2(i-t1)=(sum(sum(Ds(:,(i-t1):(i+t1))).^2)).^al;
    end
    for i=1:N
        Re0(k,i)=1/(1-al)*log2((rb1(:,i)./rb2(:,i))+eps);
    end
    
    
    Ds = abs(SST);
    Ds = [fliplr(Ds(:,1:t1)) Ds fliplr(Ds(:,end-t1+1:end))]; %padding
    for i=t1+1:N-t1
        rb1(i-t1)=sum(sum(Ds(:,(i-t1):(i+t1)).^alfa));
        rb2(i-t1)=(sum(sum(Ds(:,(i-t1):(i+t1))).^2)).^al;
    end
    for i=1:N
        Re1(k,i)=1/(1-al)*log2((rb1(:,i)./rb2(:,i))+eps);
    end
    
    Ds = abs(SST2);
    Ds = [fliplr(Ds(:,1:t1)) Ds fliplr(Ds(:,end-t1+1:end))]; %padding
    for i=t1+1:N-t1
        rb1(i-t1)=sum(sum(Ds(:,(i-t1):(i+t1)).^alfa));
        rb2(i-t1)=(sum(sum(Ds(:,(i-t1):(i+t1))).^2)).^al;
    end
    for i=1:N
        Re2(k,i)=1/(1-al)*log2((rb1(:,i)./rb2(:,i))+eps);
    end
      
end

[x,y] = min(Re0,[],1);
sig_u = sgm(y);          

[x,y] = min(Re1,[],1);
sig_sst1 = sgm(y);          

[x,y] = min(Re2,[],1);
sig_sst2 = sgm(y);          
end



function [STFT SST VSST] = stft_sst2(s,gamma,sigma)
%  stft_sst2 : computes the STFT of a signal and different versions of synchrosqueezing/reassignment.
%   Uses a Gaussian window.

n = length(s);
nv = log2(n);
if mod(nv,1)~=0
    warning('The signal is not a power of two, truncation to the next power');
    s = s(1:2^floor(nv));
end
n = length(s);
s = s(:);


% Optional parameters
if nargin<5
   bt = 1:n;
   ft = 1:n/2;
end
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
gpp = (-1/sigma^2+1/sigma^4*t.^2) .* g; % g''

% Initialization
STFT = zeros(neta,nb);
SST = zeros(neta,nb);
RSTFT = zeros(neta,nb);
OSST = zeros(neta,nb);
VSST = zeros(neta,nb);
GAR = zeros(neta,nb);
omega = zeros(neta,nb);
tau = zeros(neta,nb);
omega2 = zeros(neta,nb);
phipp = zeros(neta,nb);

%% Computes STFT and reassignment operators
df = ft(2)-ft(1);
db = bt(2)-bt(1);
for b=1:nb
	% STFT, window g
    tmp = (fft(x(bt(b):bt(b)+n-1).*g))/n;
    vg = tmp(ft);
    
    % STFT, window gp
    tmp = fft(x(bt(b):bt(b)+n-1).*gp)/n;
    vgp = tmp(ft);
    
    % operator omega
    omega(:,b) = (ft-1)'-real(tmp(ft)/2/1i/pi./vg);
    
    % STFT, window xg
    tmp = fft(x(bt(b):bt(b)+n-1).*t.*g)/n;
    vxg = n*tmp(ft);
    
    % operator tau
    tau(:,b) = n*real(tmp(ft)./vg);

    % STFT, window xxg
    tmp = fft(x(bt(b):bt(b)+n-1).*t.^2.*g)/n;
    vxxg = n^2*tmp(ft);
    
    % STFT, window gpp
    tmp = fft(x(bt(b):bt(b)+n-1).*gpp)/n;
    vgpp = tmp(ft);
    
    % STFT, window vxgp
    tmp = fft(x(bt(b):bt(b)+n-1).*t.*gp)/n;
    vxgp = n*tmp(ft);
        
    % operator hat q: estimation of frequency modulation
	phipp(:,b) = real(1/2/1i/pi*(vgpp.*vg-vgp.^2)./(vg.^2+vxg.*vgp-vxgp.*vg));
     
    % Second-order instantaneous frequency
    omega2(:,b) = omega(:,b) - phipp(:,b).*tau(:,b);
    
    % Storing STFT
    STFT(:,b) = vg.* exp(1i*pi*(ft-1)'); % compensates the tranlation 1/2 of s
        
end

% Reassignment step
for b=1:nb
    for eta=1:neta
        if abs(STFT(eta,b))>gamma
            k = 1+round((omega(eta,b)-ft(1)+1)/df);
            if k>=1 && k<=neta
                SST(k,b) = SST(k,b) + STFT(eta,b);
                l = round((tau(eta,b))/db);
            end
            k = 1+round((omega2(eta,b)-ft(1)+1)/df);
            if k>=1 && k<=neta
                VSST(k,b) = VSST(k,b) + STFT(eta,b);
            end
        end
    end
end

STFT = df*STFT;
SST = df*SST;
VSST = df*VSST;
end



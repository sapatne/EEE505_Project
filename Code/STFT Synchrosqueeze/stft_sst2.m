function [STFT SST VSST] = stft_sst2(s,gamma,sigma)
%  STFT_SST2 computes the STFT, SST and VSST of a signal 
% using a Gaussian window.
%
% INPUTS:  
%   s: real or complex signal, must be of length 2^N
%   gamma: threshold
%   sigma: window parameter
%
% INPUTS:   
%   STFT: the short-time Fourier transform
%   SST: standard synchrosqueezing
%   VSST: vertical second-order synchrosqueezing
%   omega: instantaneous frequency (vertical reassignment operator)
%   tau: phase delay (horizontal reassignment operator)
%   omega2: second-order instantaneous frequency



% checking length of signal
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
gpp = (-1/sigma^2+1/sigma^4*t.^2) .* g; % g''

% Initialization
STFT = zeros(neta,nb);
SST = zeros(neta,nb);
VSST = zeros(neta,nb);
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
    %phipp(:,b) = real(1/2/1i/pi*(vgpp.*vg-vgp.^2)./(vxg.*vgp-vxgp.*vg));
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
                % Vertical reassignment: SST 
                SST(k,b) = SST(k,b) + STFT(eta,b);
            end
            k = 1+round((omega2(eta,b)-ft(1)+1)/df);
            if k>=1 && k<=neta
                % second-order Vertical reassignment: VSST
                VSST(k,b) = VSST(k,b) + STFT(eta,b);
            end
        end
    end
end

STFT = df*STFT;
SST = df*SST;
VSST = df*VSST;

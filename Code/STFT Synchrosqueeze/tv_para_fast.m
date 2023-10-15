function vt = tv_para_fast(s,sgm_1,sgm_u,d_sgm,lmd,gamma)
% TV_PARA_FAST calculates the time-varying parameter \sigma based on the
% well-separation condition.
%
% This is a fast version, supposing we know the number of components and
% other information.
%
% Inputs:
%    s:   the signal to be analyzed
%    sgm_1,sgm_u:   the possible range of the parameter sigma
%    sgm_d:   increment
%    lmd:   a parameter to determine the length of a window
%    gamma:   a threshold
% 
%  Outputs:   
%    vt:  the time-varying parameter
%
% Example:
% sgm_1 = 0.005, d_sgm = 0.001;
% sgm = sgm_1:d_sgm:0.1;
% gamma1 = 0.3;
% [sgm_u,sgm_R1,sgm_R2] = Renyi_entropy_STFT_SST(s,sgm,al,t1);% Entropy
% lmd = 1/5;    %duration
% ci_est = tv_parameter_fast(s,sgm_1,sgm_u,d_sgm,lmd,gamma1);  % Alg.1: etimate the time-varying parameter
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
vt = sgm_u;
nc = 2;    %  number of components  
rk_rr = [64,50];% just for test;

arfa = 1/(2*pi)*sqrt(2*log(1/lmd));
n = length(s);
vt = sgm_u;

%%%% Estimate the time-varying parameter
for kk = 1:n 
    ci = sgm_u(kk);
    flag = 1;
    F1 = stft_ti(s,ci,kk);
    a_1 = local_max(F1,gamma);
    la1 = length(a_1);
    gk = zeros(1,la1);
    hk = zeros(1,la1);
    ck = zeros(1,la1);
    rk = zeros(1,la1);
    
    while la1>2 && ci > sgm_1   %% for tow components
        ci = ci-d_sgm;
        F1 = stft_ti(s,ci,kk);
        a_1 = local_max(F1,gamma);
        la1 = length(a_1);
    end
    
    if la1 == 2
    while ci > sgm_1 && flag ==1
        vt(kk) = ci;
        ck = a_1;
        gk = ck-arfa/ci;   % to save computation, we use the harmonic mode first      
        hk = ck+arfa/ci;   % because the duration of harmonic is smaller than LFM
        if sum(gk(2:end)>hk(1:end-1))~= la1-1
            flag = 0;
        else
        rk = rk_rr;% for test;
            for ll = 1:la1
                ck(ll) = a_1(ll);
                if abs(rk(ll)) <= 3
                    gk(ll) = ck(ll)-arfa/ci; 
                    hk(ll) = ck(ll)+arfa/ci; 
                else
                    gk(ll) = ck(ll)-arfa*(1/ci+2*pi*ci*rk(ll));
                    hk(ll) = ck(ll)+arfa*(1/ci+2*pi*ci*rk(ll));
                end
            end
            if sum(gk(2:end)>hk(1:end-1))~= la1-1
                flag = 0;
            else
                ci = ci-d_sgm;
                F1 = stft_ti(s,ci,kk);
                a_2 = local_max(F1,gamma);
                la2 = length(a_2);
                if la1 ~= la2
                    flag = 0;
                end
            end
        end
    end
    elseif la1 >2
    while ci > sgm_1 && flag ==1
        vt(kk) = ci;
        ak_r = a_1;
        ck = ak_r;
        gk = ck-arfa/ci;   % to save computation, we use the harmonic mode first      
        hk = ck+arfa/ci;   % because the duration of harmonic is smaller than LFM
        if sum(gk(2:end)>hk(1:end-1))~= la1-1
            flag = 0;
        else
            for ll = 1:la1
                rk(ll) = ch_ra_esti(s,a_1(ll),kk,ci,lmd);%    
                ck(ll) = a_1(ll);
                if abs(rk(ll)) <= 5e-5
                    gk(ll) = ck(ll)-arfa/ci; 
                    hk(ll) = ck(ll)+arfa/ci; 
                else
                    gk(ll) = ck(ll)-arfa*(1/ci+2*pi*ci*rk(ll));
                    hk(ll) = ck(ll)+arfa*(1/ci+2*pi*ci*rk(ll));
                end
            end
            if sum(gk(2:end)>hk(1:end-1))~= la1-1
                flag = 0;
            else
                ci = ci-d_sgm;
                F1 = stft_ti(s,ci,kk);
                a_2 = local_max(F1,gamma);
                la2 = length(a_2);
                if la1 ~= la2
                    flag = 0;
                end
            end
        end
    end
    
    end
end

end







    









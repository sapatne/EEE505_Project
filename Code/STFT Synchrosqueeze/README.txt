  
TOOLBOX of adaptive_FSST_v1

   1. This is a collection of Matlab routines developed for the adaptive STFT, adaptive FSST and 2nd-order adaptive FSST described in [1], and the entropy-based adaptive FSST with regular phase transformation in [2].


   2. List of functions included:

    MAIN ROUTINES:

    MAIN - This is an example for the two-component LFM signal

    ADAP_STFT_SST - Computes the adaptive STFT, the adaptive FSST

    STFT_SST2 - Computes STFT, FSST and VSST (2nd-order conventional FSST) of a signal using a Gaussian window
     
    REGULAR_PT_ADAP_SST- Computes the regular-phase-transform-based adaptive FSST in [2]

    REGULAR_PT_ADAP_SST2 - Computes the regular-phase-transform-based 2nd-order adaptive FSST in [2] 
      
    
    MAJOR ROUTINES:

    Renyi_STFT_SST - Finds the time-varying sigma by Renyi entropy

    IMAGESQ - Displays time-frequency result of the synchrosqueezing transform
    
    Renyi_ENTROPY_GLOBAL - Computes concentration of a time-frequency distribution.

    STFT_TI - Computes the STFT for fixed sigma and t 

    CH_RA_ESTI - Estimates the chirp rate of a component

    LOCAL_MAX - Finds the local maxima


  3. ACKNOWNLEGEMENT: We used the codes for the 2nd-order conventional FSST (called VSST) from the toolbox FSSTn-master in our STFT_SST2. We thank Duong Hung PHAM and Sylvain Meignen for their toolbox FSSTn-master. 


  4. REFERENCES
     [1] Lin Li, Haiyan Cai, Hongxia Han, Qingtang Jiang, and Hongbing Ji, "Adaptive short-time Fourier transform and synchrosqueezing transform for non-stationary signal separation," preprint, 2018, arXiv:1812.11292.

     [2] Y.-L. Sheu, L.-Y. Hsu, P.-T. Chou, and H.-T. Wu, 'Entropy-based time-varying window width selection for nonlinear-type time-frequency analysis," Int'l J. Data Sci. Anal., 3 (2017), 231-245. 

     [3] T. Oberlin, S. Meignen, and V. Perrier,"Second-order synchrosqueezing transform or invertible reassignment? Towards ideal time-frequency representations," IEEE Trans. Signal Proc., 63 (2015), 1335-1344.

     [4] Lin Li, Haiyan Cai, and Qingtang Jiang, "Adaptive synchrosqueezing transform with a time-varying parameter for non-stationary signal separation," preprint, 2018, arXiv:1812.11364.

     [5] Haiyan Cai, Qingtang Jiang, Lin Li, and Bruce W. Suter, "Analysis of adaptive short-time Fourier transform-based synchrosqueezing," preprint, 2018, arXiv:1812.11033.


  5. Copy adaptive_FSST_v1 is copyright reserved. For further information, please contact Lin Li at  lilin@xidian.edu.cn or Qingtang Jiang at jiangq@umsl.edu.

May 2018 


This toolbox can not be used for commercialization without the authorization of its author(s). If you use this toolbox for research, please must cite the following paper:
[1] Lin Li, Haiyan Cai and Qingtang Jiang, "Adaptive synchrosqueezing transform with a time-varying parameter for non-stationary signal separation," preprint, 2018, arXiv:1812.11364.
[2] Lin Li, Haiyan Cai, Hongxia Han, Qingtang Jiang, and Hongbing Ji, "Adaptive short-time Fourier transform and synchrosqueezing transform for non-stationary signal separation," preprint, 2018, arXiv:1812.11292.
[3] Haiyan Cai, Qingtang Jiang, Lin Li, and Bruce W. Suter, "Analysis of adaptive short-time Fourier transform-based synchrosqueezing," preprint, 2018, arXiv:1812.11033.


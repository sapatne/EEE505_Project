function signal = signal_gen(type, snr)
%   Generate signals for use in EEE505 project
%   Can generate 3 different signals: linear chirp, quadratic chirp,
%   multicomponent signal
%   
%   INPUTS
%   type = signal type ('l' = linear, 'q' = quadratic, 'm' =
%   multicomponent)
%   disp = boolean to show a spectrogram of the generated signal
%   
%   OUTPUT
%   signal = generated signal
%
    signal = 0;
    type = 'm';     %For purposes of running this file individually
    snr = 10;       %For purposes of running this file individually, Larger snr leads to less noise
    fs = 16384;
    Td = 1;
    N = Td*fs + 1;
    t = linspace(0, Td, N);
    F = 1:fs/2;
    switch type
        case 'l'
            f0 = 7000;
            f1 = 1000;
            signal = chirp(t, f0, Td, f1);
            signal = awgn(signal,snr,'measured');
            name = "LinearChirp";
        case 'q'
            fo = 6500;
            f1 = 2500;
            signal = chirp(t,fo,Td,f1,'quadratic',[],'convex');
            signal = awgn(signal,snr,'measured');
            name = "QuadraticChirp";
        case 'm'
            % Generate linear chirp
            f01 = 3500;
            f11 = 6750;
            x1 = chirp(t, f01, Td, f11);
    
            % Generate quadratic chirp
            f02 = 3000;
            f12 = 250;
            x2 = chirp(t, f02, Td, f12, 'quadratic', [], 'convex');

            % Combine both
            signal = x1 + x2;
            signal = awgn(signal,snr,'measured');
            name = "MultiChirp";
        otherwise
            error("ERROR: Not a valid signal type.");
    end
    
    [s1, f1, t1] = stft(signal, fs, "Window", hanning(301), "FFTLength", 1024, "FrequencyRange", "onesided");
    
    figure;
    imagesc(t1, f1, abs(s1));
    axis xy;
    plotname = "Plots/Signals/" + name + ".fig";
    savefig(plotname);

    filename = "Data/" + name + ".mat";
    save(filename, "signal");
end
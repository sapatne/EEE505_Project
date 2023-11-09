function signal = signal_gen(type, snr)
%   Generate signals for use in EEE505 project
%   Can generate 3 different signals: linear chirp, quadratic chirp,
%   multicomponent (linear + quadratic) chirp signal
%   
%   INPUTS
%   type = signal type, (DEFAULT = 'm') 
%          ('l'/'q'/'m' = linear/quadratic/multicomponent)
%   snr = desired SNR of the generated signal (scalar), (DEFAULT = 10)
%   
%   OUTPUT
%   signal = generated signal of specified type and SNR
%
    
    % Set default arguments
    arguments
        type = 'm'
        snr = 10 % Output SNR --> High SNR = less Noise
    end



    fs = 16384; % In powers of 2
    Td = 1;
    N = Td * fs + 1;
    t = linspace(0, Td, N);

    switch type
        case 'l'
            f0 = 7000; % Start
            f1 = 1000; % End
            signal = chirp(t, f0, Td, f1);
            signal = awgn(signal, snr, 'measured');
            name = "LinearChirp";
            titleplt = "Linear Chirp Generated Signal";
        case 'q'
            fo = 6500; % Start
            f1 = 2500; % End
            signal = chirp(t, fo, Td, f1, 'quadratic', [], 'convex');
            signal = awgn(signal, snr, 'measured');
            name = "QuadraticChirp";
            titleplt = "Quadratic Chirp Generated Signal";
        case 'm'
            % Generate real linear chirp
            f01 = 3500; % Start
            f11 = 6750; % End
            x1 = chirp(t, f01, Td, f11); 
    
            % Generate real quadratic chirp
            f02 = 3000; % Start
            f12 = 250; % End
            x2 = chirp(t, f02, Td, f12, 'quadratic', [], 'convex');

            % Combine both and add noise
            signal = x1 + x2;
            signal = awgn(signal, snr, 'measured');
            name = "MultiChirp";

            titleplt = "Multi-Component Generated Signal";
        otherwise
            msg = sprintf("ERROR: Not a valid signal type.\nCheck help for valid arguments.");
            error(msg);
    end
    
    % Get one-sided Spectrogram for real signal
    [s1, f1, t1] = stft(signal, fs, "Window", hanning(301), "FFTLength", 1024, "FrequencyRange", "onesided");
    
    figure;
    imagesc(t1, f1, abs(s1));
    axis xy;
    title(titleplt, FontSize=15);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    plotname = "Plots/Signals/" + name;
    savefig(plotname);
    saveas(gcf, plotname, 'png');

    filename = "Data/" + name + ".mat";
    save(filename, "signal");
end
function [s, f, r, c] = signal_gen(type, noise, N)
%   Generate signals for use in EEE505 project
%   Can generate 3 different signals: linear chirp, quadratic chirp,
%   multicomponent (linear + quadratic) chirp signal
%   
%   INPUTS
%   N = sample length/sampling frequency
%   t = signal time vector
%   type = signal type, (DEFAULT = 'n') 
%          ('n'/'o' = non-overlapping/overlapping multicomponent signal)
%   noise = set if generated signal to have AWGN, (DEFAULT = false)
%   
%   OUTPUT
%   signal = generated signal of specified type and noise
%
    
    % Set default arguments
    arguments
        type = 'n'
        noise = false 
        N = 256
    end
    
    t = linspace(0, 1-1/N, N);
    switch type
        case 'n'
            % Generate real linear chirp
            c1 = 50; % Start
            c2 = 102; % End
            x1 = chirp(t, c1, 1, c2); 
            r1 = (c2 - c1);
            f1 = c1 + r1*t;

            % Generate real quadratic chirp
            d1 = 34; % Start
            d2 = 12; % End
            % x2 = chirp(t, d1, 1, d2, 'quadratic');
            % r2 = d2 - d1
            % f2 = d1 + r2*t.^2;
            x2 = chirp(t, d1, 1, d2); 
            r2 = (d2 - d1);
            f2 = d1 + r2*t;

            % Combine both and add noise
            s = x1 + x2;
            if (noise)
                s = awgn(s, snr, 'measured');
            end
            r = [r1; r2];
            c = [c1; d1];
            f = [f1; f2];
        case 'o'
            % Generate real linear chirp
            c1 = 12; % Start
            c2 = 62; % End
            x1 = chirp(t, c1, 1, c2); 
            r1 = (c2 - c1);
            f1 = c1 + r1*t;
    
            % Generate real linear chirp
            d1 = 98; % Start
            d2 = 34; % End
            x2 = chirp(t, d1, 1, d2);
            r2 = (d2 - d1);
            f2 = c1 + r2*t;

            % Combine both and add noise
            s = x1 + x2;
            if (noise)
                s = awgn(s, snr, 'measured');
            end

            r = [r1; r2];
            c = [c1; d1];
            f = [f1; f2];

        otherwise
            msg = sprintf("ERROR: Not a valid signal type.\nCheck help for valid arguments.");
            error(msg);
    end

end
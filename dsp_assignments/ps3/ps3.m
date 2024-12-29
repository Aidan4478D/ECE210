% SPDX-License-Identifier: GPL-3.0-or-later
%
% ECE310 Digital Signal Processing - Problem Set 1
% Copyright (C) 2024 Aidan Cusa <aidancusa@gmail.com>, Kristof Jablonowski,
% Noam Schuck

% Load data and set common parameters
close all;
clc;
load projIB.mat

% Global parameters
fs = 44100;          % Sampling frequency
fs_stop = 4000;      % Stopband frequency
rp = 3;              % Passband ripple (40-37 = 3 dB)
rs = 93.5;           % Stopband attenuation in dB
gain = 10^(38.5/20); % Convert gain from dB
nyq = fs/2;          % Nyquist frequency

filter_types = {'Butterworth', 'Chebyshev1', 'Chebyshev2', 'Elliptic', 'Parks-McClellan', 'Kaiser'};

for i = 1:length(filter_types)
    fprintf('\n%s Filter Analysis:\n', filter_types{i});
    
    % Handle FIR and IIR filters differently
    if strcmp(filter_types{i}, 'Parks-McClellan') || strcmp(filter_types{i}, 'Kaiser')
        is_fir = true;
        M = 4;
        fp = 100;

        if strcmp(filter_types{i}, 'Parks-McClellan')
            Dpass = 0.10099735734;
            Dstop = 2.1134890398e-05;
            N = 16;
            dens = 20;
            [n, Fo, Ao, W] = firpmord([fp, fs_stop]/(nyq), [1 0], [Dpass, Dstop]);
            b = firpm(n, Fo, Ao, W, {dens});
            a = 1;
            b = b * (gain+50);
        else % Kaiser
            Dpass = 0.17099735734;
            Dstop = 2.1134890398e-05;
            N = 16;
            [n, wn, beta, type] = kaiserord([fp fs_stop]/(nyq), [1 0], [Dstop Dpass]);
            b = fir1(n, wn, type, kaiser(n+1, beta));
            a = 1;
            b = b * gain;
        end
    else
        % for IIR filters find optimal passband frequency
        is_fir = false;
        [fp, b, a] = find_optimal_fp(filter_types{i}, fs_stop, nyq, rp, rs, gain);
    end
    

    [h, w] = freqz(b, a, 1024);
    figure;
    subplot(2,1,1);
    plot(w/pi * nyq, 20*log10(abs(h)));
    grid on;
    title(sprintf('%s Filter Response', filter_types{i}));
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    ylim([-100 45]);
    
    if is_fir
        M = 4; % Downsampling factor
        % polyphase filtering with downsampling using this sexy function
        downsampled_signal = upfirdn(noisy, b, 1, M);
    else
        filtered_signal = filter(b, a, noisy);
        % Downsample by 4 (compressor)
        downsampled_signal = filtered_signal(1:4:end);
    end
    
    mults = count_multiplies(b, a, is_fir);
    fprintf('Passband frequency: %.1f Hz\n', fp);
    fprintf('Number of multiplies: %.1f\n', mults);
    

    subplot(2,1,2);
    plot((1:5000)/1000, noisy(1:5000), 'b', ...
         (1:length(downsampled_signal))/ (fs/4), downsampled_signal, 'r');
    grid on;
    title('Time Domain Response');
    xlabel('Time (s)');
    ylabel('Amplitude');
    legend('Original', 'Filtered');
    
    % Play audio
    sound(downsampled_signal, fs/4);
    pause(6);
    clear sound;
    pause(1);
end

% https://www.mathworks.com/help/signal/ref/upfirdn.html
% Total Multiplications = (Lh - 1) * Lx / q
% Multiplications per Input Sample = (Lh - 1) / q​
    % h(n) corresponds to the FIR filter coefficients b
function mult_count = count_multiplies(b, a, is_fir)
    if is_fir
        % FIR with polyphase is valid
        M = 4; % Downsampling factor
        mult_count = ceil((length(b) - 1) / M);
    else
        % each coefficient in a (excluding a(1), which is typically normalized to 1) requires one multiplication per output sample.
        mult_count = length(b) + length(a) - 1;  % -1 because a[0] is always 1
    end
end

% Function to find optimal passband frequency
function [opt_fp, opt_b, opt_a] = find_optimal_fp(filter_type, fs_stop, nyq, rp, rs, gain)
    fp_min = 1000;  % Minimum passband frequency
    fp_max = fs_stop;
    opt_fp = fp_min;
    opt_b = [];
    opt_a = [];
    
    while (fp_max - fp_min) > 1
        fp_test = (fp_min + fp_max)/2;
        wp_test = fp_test/nyq;
        ws = fs_stop/nyq;
        
        switch filter_type
            case 'Butterworth'
                [n, wn] = buttord(wp_test, ws-500/nyq, rp, 55);
                [b, a] = butter(n, wn);
                is_fir = false;
            case 'Chebyshev1'
                [n, wn] = cheb1ord(wp_test, ws, rp, rs);
                [b, a] = cheby1(n, rp, wn);
                is_fir = false;
            case 'Chebyshev2'
                [n, wn] = cheb2ord(wp_test, ws, rp, rs);
                [b, a] = cheby2(n, rs, wn);
                is_fir = false;
            case 'Elliptic'
                [n, wn] = ellipord(wp_test, ws, rp, rs);
                [b, a] = ellip(n, rp, rs, wn);
                is_fir = false;
        end
        
        % Apply gain
        b = b * gain;
        
        % Check if multiplication constraint is met
        if count_multiplies(b, a, is_fir) <= 17
            fp_min = fp_test;
            if fp_test > opt_fp
                opt_fp = fp_test;
                opt_b = b;
                opt_a = a;
            end
        else
            fp_max = fp_test;
        end
    end
end

% SPDX-License-Identifier: GPL-3.0-or-later
%
% ECE310 Digital Signal Processing - Problem Set 1
% Copyright (C) 2024 Aidan Cusa <aidancusa@gmail.com>, Kristof Jablonowski,
% Noam Schuck

[x, Fs_original] = audioread('wagner.wav');

% Perform sampling rate conversion
y = srconvert(x);

% Write the output file
audiowrite('wagner_24k.wav', y, 24000);

% Display completion message
disp('Conversion complete. Output saved as wagner_24k.wav');

% Test with a zeros input
y2 = srconvert([1 zeros(1, 3000)]);
figure(11)
verify(y2, "Zeros")

function y = upsample_polyphase(x, E, phase_shift)
    L = size(E, 1);  % Upsampling factor
    y = zeros(1, L * length(x));  % Allocate space for upsampled signal
    phase_shift = mod(phase_shift, L);  % Ensure the phase shift is within the range [0, L-1]
    
    for i = 1:L
        shifted_idx = mod(i + phase_shift - 1, L) + 1;  % Apply phase shift
        y(shifted_idx:L:end) = fftfilt(E(i,:), x);  % Use FFT-based convolution for filtering
    end
end

function y = downsample_polyphase(x, M, phase_shift)
    % Downsampling without filtering, only using polyphase components
    num_samples = ceil((length(x) - phase_shift) / M);  % Calculate the number of output samples
    y = zeros(1, num_samples);  % Allocate space for downsampled signal
    phase_shift = mod(phase_shift, M);  % Ensure the phase shift is within the range [0, M-1]
    
    idx = 1;  % Initialize index for output y
    for i = phase_shift+1:M:length(x)
        if idx <= length(y)
            y(idx) = x(i);  % Downsample by selecting polyphase components
            idx = idx + 1;
        else
            break;
        end
    end
end


% Local function below

function y = srconvert(x)
    % 5x upsampling with a phase shift of 0 (no shift)
    h_5x = design_filter(0.2, 0.24, 0.001, 1);  % Design a filter for upsampling
    E_5x = poly1(h_5x, 5);  % Create polyphase components for 5x upsampling
    y = upsample_polyphase(x, E_5x, 0);  % No phase shift for the first stage
    
    % Verify the 5x upsampling filter
    figure(1);
    verify(h_5x, '5x Upsampling Filter');
    
    % Six 2x upsampling stages with progressive phase shifts
    for i = 1:6
        h_2x = design_filter(0.5, 0.6, 0.001, 1);  % Design filter for 2x upsampling
        E_2x = poly1(h_2x, 2);  % Create polyphase components for 2x upsampling
        phase_shift = mod(i, 2);  % Apply a small phase shift that alternates between 0 and 1
        y = upsample_polyphase(y, E_2x, phase_shift);
        
        % Verify each 2x upsampling filter
        figure(i + 1);
        %verify(h_2x, ['2x Upsampling Filter Stage ' num2str(i)]);
    end
    
    % Two 7x downsampling stages with phase shifts (no filtering)
    for i = 1:2
        phase_shift = mod(i, 7);  % Apply small phase shifts
        y = downsample_polyphase(y, 7, phase_shift);  % No filtering, just polyphase selection
        
        % Verify each 7x downsampling stage
        figure(7 + i);
        %verify(y, ['7x Downsampling Stage ' num2str(i)]);
    end
    
    % Final 3x downsampling with a phase shift of 1 (no filtering)
    y = downsample_polyphase(y, 3, 1);  % Apply a small phase shift
    
    % Verify the 3x downsampling stage
    figure(10);
    verify(y, '3x Downsampling Stage');
    
    % Adjust gain if needed
    y = y * 150;
end


function h = design_filter(Fpass, Fstop, Dpass, Dstop)
    % Calculate the order from the parameters using FIRPMORD
    [N, Fo, Ao, W] = firpmord([Fpass, Fstop], [1 0], [Dpass, Dstop]);
    
    % Calculate the coefficients using the FIRPM function
    % Note: We're using a default density factor here. Adjust if needed.
    dens = 20;
    b = firpm(N, Fo, Ao, W, {dens});
    
    % Create a discrete-time FIR filter object
    Hd = dfilt.dffir(b);
    
    % Extract the filter coefficients
    h = Hd.Numerator;
end

function yes = verify(h, title_str)
    [R, G, A] = examlpf(h, 147/320, 147/320*1.2);
    
    % Passband magnitude
    yes = R <= 0.1;
    subplot(3,1,1);
    if R <= 0.1
        title([title_str ': Passband Response, Ripple = ' num2str(R) ' dB, OK']);
    else
        title([title_str ': Passband Response, Ripple = ' num2str(R) ' dB, too big']);
    end
    
    % Group delay
    subplot(3,1,2);
    if G <= 720
        title([title_str ': Group Delay in Passband, max-min = ' num2str(G) ', OK']);
    else
        title([title_str ': Group Delay in Passband, max-min = ' num2str(G) ', too big']);
    end
    yes = yes * (G <= 720);
    
    % Stopband magnitude
    subplot(3,1,3);
    if A <= -70
        title([title_str ': Overall Response, Attenuation = ' num2str(A) ' dB, OK']);
    else
        title([title_str ': Overall Response, Attenuation = ' num2str(A) ' dB, too small']);
    end
    yes = yes * (A <= -70);
    
    fprintf('%s:\n', title_str);
    fprintf('Passband Ripple: %5.3f dB \n', R);
    fprintf('Group delay Variation: %5d samples \n', G);
    fprintf('Stopband Attenuation: %5.3f dB \n\n', A);
end
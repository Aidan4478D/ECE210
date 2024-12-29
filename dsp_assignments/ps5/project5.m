% Aidan Cusa, Kristof Jablonowski, Noam Schuck

clear; 
close all; 
clc;

[x, fs] = audioread('bjruf.wav');

N = 1024; % FFT size
L = N; % Window length = FFT size

%R = L/2; % uncomment if want to use barlett window (comment if rect)
R = L; % uncomment if want to use rect window (comment if barlett)

numBlocks = ceil(length(x)/R);
x = x(:)'; % make x a row vector

% For no overlap, we did:
% x = [x zeros(1, numBlocks*R - length(x))];

% need at least (numBlocks-1)*R + L samples total to ensure we have enough padding
requiredLength = (numBlocks-1)*R + L;
x = [x zeros(1, requiredLength - length(x))];

w = rectwin(L)'; % uncomment if want to use rect window (comment if barlett)
%w = bartlett(L)'; % uncomment if want to use barlett window (comment if rect)

X_r = zeros(numBlocks, N);  
speech_presence = zeros(numBlocks, 1);

for r = 1:numBlocks
    start_idx = (r-1)*R + 1;
    seg = x(start_idx:start_idx+L-1) .* w; 
    X_block = fft(seg, N); % forward DFT
    X_r(r,:) = X_block;
end

% Speech Detection Heuristic - measure energy in frequency domain and threshold it.
% if the energy is above a certain fraction of the peak energy, it assumes that block contains speech.

% block compute total energy
block_energies = sum(abs(X_r).^2, 2);

% highest-energy block is taken as reference, and any block whose energy is greater than 10% of the maximum observed energy is considered to contain speech
energy_threshold = 0.1 * max(block_energies);

% 0 if no speech, 1 if speech
speech_presence(block_energies > energy_threshold) = 1;

x_hat = zeros(1, (numBlocks-1)*R + L);
for r = 1:numBlocks
    % IDFT (implements 10.43)
    x_r_block = ifft(X_r(r,:), N); 
    x_r_block = real(x_r_block(1:L));  % apply window w[m] to recovered segment

    % overlap add (implements 10.44)
    start_idx = (r-1)*R + 1;
    x_hat(start_idx:start_idx+L-1) = x_hat(start_idx:start_idx+L-1) + x_r_block;
end

% trim to original length
x_hat = x_hat(1:length(x));

% verify reconstruction quality
reconstruction_error = norm(x_hat - x) / norm(x);
disp(['Reconstruction error: ' num2str(reconstruction_error)]);

% create a 0-1 vector aligned with the original signal for visualization
speech_binary = zeros(size(x));
for r = 1:numBlocks
    start_idx = (r-1)*R + 1;
    end_idx = start_idx + L - 1;
    if end_idx > length(speech_binary)
        end_idx = length(speech_binary);
    end
    speech_binary(start_idx:end_idx) = speech_presence(r);
end

t = (0:length(x)-1)/fs;
figure;
subplot(2,1,1);
plot(t, x);
xlabel('Time (s)');
ylabel('Amplitude');
title('Original Audio Signal');

subplot(2,1,2);
plot(t, speech_binary, 'r');
xlabel('Time (s)');
ylabel('Speech Activity');
title('Detected Speech Activity (1=Speech, 0=No Speech)');
ylim([-0.1 1.1]);

% play the reconstructed sound
soundsc(x_hat, fs);
pause(10);
clear sound;

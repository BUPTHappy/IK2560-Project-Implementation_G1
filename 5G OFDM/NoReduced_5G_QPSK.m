% Parameters for mixed numerology and different NFFT values
Nfft_values = [128, 256, 512]; % Different NFFT sizes for comparison
subcarriers = [72, 180, 300]; % Different subcarrier numerologies (e.g., 15 kHz, 30 kHz, 60 kHz)
SNR = 20; % Signal-to-noise ratio in dB
num_symbols = 1000; % Number of OFDM symbols to simulate
mod_order = 4; % QPSK modulation order (QPSK has 4 symbols)

% Generate random data for each numerology
data = cell(1, length(subcarriers));
modulated_data = cell(1, length(subcarriers));
transmitted_signal = cell(1, length(subcarriers));
received_signal = cell(1, length(subcarriers));

% Generate data for each numerology
for i = 1:length(subcarriers)
    data{i} = randi([0 mod_order-1], num_symbols, subcarriers(i)); % Random data symbols
    modulated_data{i} = qammod(data{i}, mod_order, 'UnitAveragePower', true); % QPSK modulation
end

% PAPR Calculation for Different NFFT
papr_values = cell(length(Nfft_values), length(subcarriers));

% Loop through different NFFT sizes
for nfft_idx = 1:length(Nfft_values)
    Nfft = Nfft_values(nfft_idx);
    
    % OFDM Modulation (Inverse FFT) for each numerology and NFFT size
    for i = 1:length(subcarriers)
        % Perform IFFT to get the time-domain signal for the current NFFT
        transmitted_signal{i} = ifft(modulated_data{i}, Nfft, 2); % Nfft is the FFT size
    end
    
    % Add cyclic prefix (CP)
    cp_length = Nfft / 8; % CP length (1/8 of Nfft is typical)
    for i = 1:length(subcarriers)
        transmitted_signal{i} = [transmitted_signal{i}(:,end-cp_length+1:end), transmitted_signal{i}];
    end
    
    % Channel (AWGN) and reception
    for i = 1:length(subcarriers)
        % Adding AWGN noise
        noise = (1/sqrt(2)) * (randn(size(transmitted_signal{i})) + 1i*randn(size(transmitted_signal{i})));
        received_signal{i} = transmitted_signal{i} + 10^(-SNR/20) * noise; % Adding noise
    end
    
    % PAPR Calculation for Each Numerology and NFFT size
    for i = 1:length(subcarriers)
        % PAPR Calculation
        PAPR = 10*log10(max(abs(transmitted_signal{i}).^2, [], 2) ./ mean(abs(transmitted_signal{i}).^2, 2));
        papr_values{nfft_idx, i} = PAPR;
    end
end

% Plot CCDF for different NFFT values with customized colors and line styles
figure;
hold on;

% Define colors for different NFFT values
colors = ['k', 'g', 'r']; % Blue for Nfft = 128, Green for Nfft = 256, Red for Nfft = 512

% Define line styles for different numerologies (subcarrier values)
line_styles = {'-', '--', ':'}; % Solid, dashed, and dotted lines for different numerologies

% Define line width for better visibility
line_width = 2; % You can increase this value for even thicker lines

% Loop through NFFT values and plot the CCDF of PAPR
for nfft_idx = 1:length(Nfft_values)
    for i = 1:length(subcarriers)
        % Sort PAPR values in ascending order
        sorted_papr = sort(papr_values{nfft_idx, i});
        
        % Calculate CCDF: the probability of PAPR exceeding a certain threshold
        ccdf = (1:length(sorted_papr)) / length(sorted_papr);
        
        % Plot CCDF for each numerology and NFFT size
        plot(sorted_papr, 1-ccdf, 'Color', colors(nfft_idx), 'LineStyle', line_styles{i}, 'LineWidth',line_width,...
            'DisplayName', sprintf('NFFT = %d, Numerology %d kHz', Nfft_values(nfft_idx), subcarriers(i)));
    end
end

xlabel('PAPR [dB]');
ylabel('CCDF (Pr[PAPR>PAPR0])');
title('PAPR CCDF for 5G Mixed Numerology OFDM System with Different NFFT Sizes(QPSK)');
legend;
grid on;
hold off;

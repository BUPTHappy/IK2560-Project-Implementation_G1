clc;
clear;
close all;

% 5G OFDM参数
subCarrierSpacing = 30e3;  % 子载波间隔 30 kHz
Lfft = 512;                % FFT大小
numSymbols = 100;          % OFDM符号数量
numSubcarriers = Lfft;     % 子载波数量
cpLength = 128;            % 循环前缀长度

% 调制方式 (例如 QPSK)
modOrder = 4;              % QPSK调制
MAP_qpsk = [1 -1 j -j];    % QPSK调制集
Phases = [1 -1];           % 相位集合 W=2

% PTS参数
V = 4;                     % 子块数量
Lchos = 16;                % 相位组合数 (W^(V-1))
Pchos = combvec(Phases,Phases,Phases,Phases).'; % 所有相位组合

% 生成随机符号
data = randi([0 modOrder-1], numSymbols, numSubcarriers);

% 将符号映射到QPSK
modulatedData = zeros(numSymbols, numSubcarriers);
for i = 1:numSymbols
    modulatedData(i, :) = MAP_qpsk(data(i, :) + 1);  % QPSK映射
end

% 添加循环前缀 (CP)
ofdmSignal = zeros(numSymbols, numSubcarriers + cpLength);
papr_original = zeros(1, numSymbols);  % 原始PAPR
papr_pts = zeros(1, numSymbols);       % PTS算法后PAPR

for i = 1:numSymbols
    % ----------------- 原始OFDM信号 ----------------- 
    x = ifft(modulatedData(i, :), Lfft);
    x_cp = [x(end-cpLength+1:end) x];
    ofdmSignal(i, :) = x_cp;
    % 计算原始PAPR
    peakPower = max(abs(x_cp).^2);
    avgPower = mean(abs(x_cp).^2);
    papr_original(i) = 10 * log10(peakPower / avgPower);
    
    % ----------------- PTS处理 ----------------- 
    A = reshape(modulatedData(i, :), V, Lfft/V);  % 分成V个子块
    a_ifft = ifft(A, [], 2);  % 对每个子块IFFT
    
    min_papr = Inf;
    for n = 1:Lchos
        phase_comb = Pchos(n, :).'; % 当前相位组合
        phase_matrix = repmat(phase_comb, 1, Lfft/V);
        x_pts = sum(a_ifft .* phase_matrix, 1);
        % 添加CP
        x_cp_pts = [x_pts(end-cpLength+1:end) x_pts];
        % 计算PAPR
        peakPower = max(abs(x_cp_pts).^2);
        avgPower = mean(abs(x_cp_pts).^2);
        papr_temp = 10 * log10(peakPower / avgPower);
        if papr_temp < min_papr
            min_papr = papr_temp;
        end
    end
    papr_pts(i) = min_papr;
end

% ----------------- 绘制CCDF图 -----------------
figure;
% 原始PAPR的CCDF
[cdf_orig, papr_vals_orig] = ecdf(papr_original);
semilogy(papr_vals_orig, 1-cdf_orig, 'b', 'LineWidth', 2); hold on;

% PTS后PAPR的CCDF
[cdf_pts, papr_vals_pts] = ecdf(papr_pts);
semilogy(papr_vals_pts, 1-cdf_pts, 'r', 'LineWidth', 2);

legend('Original OFDM', 'PTS Processed OFDM');
title('PAPR CCDF Comparison for 5G OFDM');
xlabel('PAPR (dB)');
ylabel('CCDF');
grid on;

% ----------------- 绘制OFDM时域信号 -----------------
figure;
plot(real(ofdmSignal(1,:)));
title('Original 5G OFDM Time Domain Signal');
xlabel('Time (samples)');
ylabel('Amplitude');
grid on;
    
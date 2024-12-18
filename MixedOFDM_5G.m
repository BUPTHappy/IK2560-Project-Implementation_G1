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

% 生成随机符号
data = randi([0 modOrder-1], numSymbols, numSubcarriers);

% 将符号映射到QPSK
modulatedData = zeros(numSymbols, numSubcarriers);
for i = 1:numSymbols
    modulatedData(i, :) = MAP_qpsk(data(i, :) + 1);  % QPSK映射
end

% 添加循环前缀 (CP)
ofdmSignal = zeros(numSymbols, numSubcarriers + cpLength);
for i = 1:numSymbols
    % IFFT将符号从频域转换到时域
    x = ifft(modulatedData(i, :), Lfft);
    
    % 循环前缀
    ofdmSignal(i, :) = [x(end-cpLength+1:end) x];
end

% 计算PAPR
PAPR = zeros(1, numSymbols);  % 存储每个符号的PAPR
for i = 1:numSymbols
    % 计算OFDM符号的时域信号
    x_t = ofdmSignal(i, :);
    
    % PAPR公式
    peakPower = max(abs(x_t).^2);  % 峰值功率
    avgPower = mean(abs(x_t).^2);  % 平均功率
    PAPR(i) = 10 * log10(peakPower / avgPower);  % 计算PAPR (dB)
end

% 绘制PAPR的CCDF图
figure;
% 绘制CCDF
[empiricalCDF, xPAPR] = ecdf(PAPR);
semilogy(xPAPR, 1-empiricalCDF, 'LineWidth', 2);
title('PAPR CCDF for 5G OFDM');
xlabel('PAPR (dB)');
ylabel('CCDF');
grid on;

% 绘制OFDM信号的时域图
figure;
plot(real(ofdmSignal(1,:)));
title('5G OFDM Time Domain Signal');
xlabel('Time (samples)');
ylabel('Amplitude');
grid on;

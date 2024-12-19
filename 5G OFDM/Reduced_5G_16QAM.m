clc;
clear;
close all;

% 5G OFDM参数
Lfft = 512;                % FFT大小
numSymbols = 1000;         % OFDM符号数量
cpLength = 128;            % 循环前缀长度
modOrder = 16;             % 16QAM调制
M = 16;                     % 调制阶数
k = log2(M);               % 每个符号的比特数

% 16QAM调制映射 (使用Gray码映射以获得更好的误码率性能)
MAP_16qam = qammod(0:M-1, M, 'gray');

Phases = [1 -1];           % 相位集合 W=2
V = 4;                     % 子块数量
Lchos = 16;                % 相位组合数 (W^(V-1))
Pchos = combvec(Phases,Phases,Phases,Phases).'; % 所有相位组合

% 生成随机比特流
data_bits = randi([0 1], numSymbols, Lfft * k);

% 将比特流分组为符号
data_symbols = reshape(data_bits, numSymbols, Lfft, k);

% 将符号映射到16QAM
modulatedData = zeros(numSymbols, Lfft);
for i = 1:numSymbols
    for j = 1:Lfft
        % 将k个比特转换为十进制数
        decimal_value = bi2de(squeeze(data_symbols(i, j, :))'); % 关键修改
        modulatedData(i, j) = MAP_16qam(decimal_value + 1);
    end
end


% 子载波间隔
subCarrierSpacings = [15e3, 30e3, 60e3];  % 子载波间隔 (15 kHz, 30 kHz, 60 kHz)
papr_original_all = zeros(length(subCarrierSpacings), numSymbols);  % 原始PAPR
papr_pts_all = zeros(length(subCarrierSpacings), numSymbols);       % PTS算法后PAPR

% 创建图形窗口
figure;

% 存储不同子载波间隔下的PAPR数据
for idx = 1:length(subCarrierSpacings)
    subCarrierSpacing = subCarrierSpacings(idx);  % 当前子载波间隔
    % 计算每个符号的PAPR
    papr_original = zeros(1, numSymbols);  % 原始PAPR
    papr_pts = zeros(1, numSymbols);       % PTS后PAPR
    for i = 1:numSymbols
        % ----------------- 原始OFDM信号 ----------------- 
        x = ifft(modulatedData(i, :), Lfft);
        x_cp = [x(end-cpLength+1:end) x];

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
    % 存储结果
    papr_original_all(idx, :) = papr_original;
    papr_pts_all(idx, :) = papr_pts;
    % 绘制每个子载波间隔下的PAPR曲线
    [cdf_orig, papr_vals_orig] = ecdf(papr_original);
    semilogy(papr_vals_orig, 1-cdf_orig,'b', 'LineWidth', 2); hold on;
    [cdf_pts, papr_vals_pts] = ecdf(papr_pts);
    semilogy(papr_vals_pts, 1-cdf_pts,'m', 'LineStyle','--', 'LineWidth', 2);
end

% 设置图形标题和标签
ylim([0.1 1]);  % 设置Y轴范围
legend('Original OFDM (mixed Numerology)','PTS Processed OFDM (mixed Numerology)');
title('PAPR CCDF Comparison for 5G OFDM with 15/30/60kH(16QAM)'); % 修改标题
xlabel('PAPR (dB)');
ylabel('CCDF');
grid on;
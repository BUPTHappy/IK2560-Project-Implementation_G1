function [finalSignal, bestPAPR] = pts_ofdm(modulatedData, Nfft, numSubblocks)
    % PTS方法减少OFDM信号的PAPR
    % 输入：
    % modulatedData: 经过QAM调制的OFDM符号，大小为 (Nfft, numSymbols)
    % Nfft: FFT大小
    % numSubblocks: 子序列的数量
    % 输出：
    % finalSignal: 经过PTS降PAPR后的OFDM信号
    % bestPAPR: 最佳PAPR值

    % 计算每个子块的大小
    subblockSize = Nfft / numSubblocks;

    % 确保Nfft能被numSubblocks整除
    if mod(Nfft, numSubblocks) ~= 0
        error('Nfft must be divisible by numSubblocks');
    end

    % 获取符号的数量
    numSymbols = size(modulatedData, 2);  % 符号数量

    % 初始化旋转因子和PAPR值
    bestPAPR = Inf;
    bestPhase = zeros(numSubblocks, 1);  % 最佳相位旋转因子
    finalSignal = zeros(Nfft, numSymbols);  % 初始化最终信号

    % 1. 遍历每个符号
    for ii = 1:numSymbols
        % 当前符号
        currentSymbol = modulatedData(:, ii);

        % 2. 将当前符号分成多个子块
        subblocks = zeros(subblockSize, numSubblocks);  % 初始化子块矩阵
        for i = 1:numSubblocks
            subblocks(:, i) = currentSymbol((i-1)*subblockSize+1:i*subblockSize);  % 按照索引切分
        end

        % 3. 遍历不同的相位旋转因子
        for phaseIdx = 1:360  % 遍历0到360度的相位
            phaseShift = exp(1j * deg2rad(phaseIdx));  % 计算当前相位旋转因子
            
            % 应用相位旋转因子到每个子块
            rotatedSubblocks = subblocks;  % 恢复子块
            for blockIdx = 1:numSubblocks
                rotatedSubblocks(:, blockIdx) = subblocks(:, blockIdx) * phaseShift;
            end
            
            % 合成信号（按列合并子块）
            combinedSignal = sum(rotatedSubblocks, 2);  % 合成信号

            % 计算PAPR
            PAPR = 10 * log10(max(abs(combinedSignal).^2) / mean(abs(combinedSignal).^2));

            % 如果当前PAPR更小，则更新最佳相位
            if PAPR < bestPAPR
                bestPAPR = PAPR;
                bestPhase = phaseShift * ones(numSubblocks, 1);  % 记录最优的相位旋转因子
            end
        end

        % 使用最佳相位旋转因子构建最终信号
        % 使用元素乘法来应用每个子块的相位因子
        for blockIdx = 1:numSubblocks
            subblocks(:, blockIdx) = subblocks(:, blockIdx) * bestPhase(blockIdx);
        end

        % 合成最终信号
        finalSignal(:, ii) = sum(subblocks, 2);  % 合成最终信号
    end
end

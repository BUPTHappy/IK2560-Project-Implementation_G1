clc;                      % 清除命令窗口
clear;                    % 清除所有变量
close all;                % 关闭所有图形窗口
warning off;              % 关闭警告信息
addpath(genpath(pwd));    % 将当前目录添加到搜索路径
rng('default');           % 设置随机数生成器的初始状态为默认值
rng(1);  

% FFT大小数组
Lfft_values = [64, 128, 256, 512, 1024];
V = 4;            % 选择的数量
MAP_qpsk  = [1 -1 j -j];  % QPSK调制集
Phases    = [1 -1];       % 相位集合 W = 2

Nbits     = 10000;      % 最大符号数量

% 产生的序列个数 = W ^(V-1) = 8
Pchos     = [1 1 1 1; 1 1 1 2; 1 1 2 1; 1 2 1 1; 2 1 1 1;...
             1 1 2 2; 1 2 1 2; 1 2 2 1; 2 2 1 1; 2 1 2 1; 2 1 1 2;...
             2 2 2 1; 2 2 1 2; 2 1 2 2; 1 2 2 2; 2 2 2 2];
Lchos     = 16;         % 选择序列的长度

% 设置图形
figure;
hold on;

% 预定义颜色
colors = lines(length(Lfft_values));  % 使用MATLAB的'lines'调色板

% 存储图例字符串
legend_entries = {};

% 循环遍历不同的FFT大小
for idx = 1:length(Lfft_values)
    Lfft = Lfft_values(idx);  % 当前FFT大小
    papr0     = zeros(1,Nbits); % 存储原始PAPR值
    papr_pts  = zeros(1,Nbits); % 存储PTS处理后的PAPR值

    for jj=1:Nbits
        Index        = floor(length(MAP_qpsk)*rand(1,Lfft))+1; % 随机选择QPSK符号索引
        X            = MAP_qpsk(Index(1,:)); % 原始频域信号
        x            = ifft(X,[],2); % 时域信号
        Pow1         = abs(x.^2); % 信号功率
        Pow2         = max(Pow1,[],2); % 峰值功率
        Pow3         = mean(Pow1,2);   % 平均功率
        papr0(jj)    = 10*log10(Pow2./Pow3); % 计算原始PAPR

        % PTS分割
        A = zeros(V,Lfft); % 初始化A矩阵
        for v=1:V
            A(v,(1+(v-1)*Lfft/4):1:v*Lfft/4) = X((1+(v-1)*Lfft/4):1:v*Lfft/4); % 将X分成V个部分放入A中
        end
        a    = ifft(A,[],2); % 对每个部分做IFFT
        % 门限函数
        Vmin = 10; % 初始化最小值
        % 寻找最优辅助信息
        for n=1:Lchos
            temp_phase = Phases(Pchos(n,:)).'; % 获取相位组合
            b = repmat(temp_phase,1,Lfft); % 复制相位组合
            c = a.*b; % 相乘
            d = sum(c); % 求和
            e = abs(d); % 取绝对值
            temp_max = max(e); % 计算最大值
            if temp_max<Vmin
                Vmin = temp_max; % 更新最小值
                Best_n = n; % 记录最佳选择
            end
        end
        aa = sum(a.*repmat(Phases(Pchos(Best_n,:)).',1,Lfft)); % 最佳选择相位调整后求和

        Pow1 = abs(aa.^2); % 新的信号功率
        Pow2   = max(Pow1,[],2); % 新的峰值功率
        Pow3   = mean(Pow1,2);   % 新的平均功率
        papr_pts(jj) = 10*log10(Pow2./Pow3); % 计算新的PAPR
    end

    % 计算累积分布函数
    [cdf1, PAPR1] = ecdf(papr0); % 原始PAPR的CDF
    [cdf2, PAPR2] = ecdf(papr_pts);     % PTS处理后的PAPR CDF

    % 绘制不同FFT大小的结果
    semilogy(PAPR1(1:end), 1-cdf1(1:end), 'Color', colors(idx,:), 'LineWidth', 2);
    semilogy(PAPR2(1:end), 1-cdf2(1:end), '--', 'Color', colors(idx,:), 'LineWidth', 2);  

    % 添加图例条目
    legend_entries{end+1} = sprintf('FFT=%d (Original)', Lfft);
    legend_entries{end+1} = sprintf('FFT=%d (PTS)', Lfft);
    
end

% 设置图例
legend(legend_entries, 'Location', 'Best');
title('Comparison of PAPR for different FFT sizes with and without PTS');
xlabel('PAPR [dB]');
ylabel('CCDF (Pr[PAPR>PAPR0])');
grid on;
hold off;

% 绘制不同FFT大小的PAPR直方图
figure;
hold on;

% 设置透明度
alpha_value = 0.5; % 透明度设置为50%

% 颜色：黄色和蓝色
color_original = [1, 1, 0];  % 黄色
color_pts = [0, 0, 1];       % 蓝色

% 使用不同颜色区分原始PAPR和PTS处理后的PAPR
for idx = 1:length(Lfft_values)
    Lfft = Lfft_values(idx);
    
    % 原始PAPR的直方图（黄色）
    histogram(papr0, 'Normalization', 'pdf', 'EdgeColor', 'none', ...
        'FaceColor', color_original, 'BinWidth', 0.1, 'FaceAlpha', alpha_value);
    
    % PTS处理后的PAPR的直方图（蓝色）
    histogram(papr_pts, 'Normalization', 'pdf', 'EdgeColor', 'none', ...
        'FaceColor', color_pts, 'BinWidth', 0.1, 'LineStyle', '--', 'FaceAlpha', alpha_value);
end

% 添加手动设置的图例
legend_entries = {'Original PAPR (Yellow)', 'PTS PAPR (Blue)'};

% 设置图例
legend(legend_entries, 'Location', 'Best');
title('PAPR Histogram for Different FFT Sizes with and without PTS');
xlabel('PAPR [dB]');
ylabel('Probability Density');
grid on;
hold off;


% 设置图形窗口
figure;
hold on;
% 使用subplot绘制多个子图
for idx = 1:length(Lfft_values)
    Lfft = Lfft_values(idx);  % 当前FFT大小
    
    % 随机选择QPSK符号索引
    Index = floor(length(MAP_qpsk)*rand(1, Lfft)) + 1;
    X = MAP_qpsk(Index(1,:)); % 频域信号
    
    % 通过IFFT计算时域信号
    x = ifft(X, Lfft);
    
    % 创建子图
    subplot(2, 3, idx);  % 2行3列的子图，第idx个子图
    plot(real(x), 'b-', 'LineWidth', 2);
    hold on;
    plot(imag(x), 'r--', 'LineWidth', 2);
    title(['OFDM Signal Time-domain (FFT=' num2str(Lfft) ')']);
    xlabel('Sample Index');
    ylabel('Amplitude');
    legend({'Real Part', 'Imaginary Part'}, 'Location', 'Best');
    grid on;
end

% 设置整体图形
sgtitle('OFDM Signal Time-domain for Different FFT Sizes');
hold off;
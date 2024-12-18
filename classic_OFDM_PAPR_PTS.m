clc;                      % 清除命令窗口
clear;                    % 清除所有变量
close all;                % 关闭所有图形窗口
warning off;              % 关闭警告信息
addpath(genpath(pwd));    % 将当前目录添加到搜索路径
rng('default');           % 设置随机数生成器的初始状态为默认值
rng(1);  


Lfft      = 128;          % FFT的大小
V         = 4;            % 选择的数量
MAP_qpsk  = [1 -1 j -j];  % QPSK调制集
Phases    = [1 -1];       % 相位集合 W = 2

Nbits     = 10000;      % 最大符号数量

% 产生的序列个数 = W ^(V-1) = 8
% 针对QPSK共有16种相位变化组合
Pchos     = [1 1 1 1; 1 1 1 2; 1 1 2 1; 1 2 1 1; 2 1 1 1;...
             1 1 2 2; 1 2 1 2; 1 2 2 1; 2 2 1 1; 2 1 2 1; 2 1 1 2;...
             2 2 2 1; 2 2 1 2; 2 1 2 2; 1 2 2 2; 2 2 2 2];
Lchos     = 16;         % 选择序列的长度


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
    % 门限函数（此处省略）
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

% 绘图
semilogy(PAPR1(1:1:end),1-cdf1(1:1:end),'-b',... % 绘制原始PAPR的CCDF
         PAPR2(1:1:end),1-cdf2(1:1:end),'-r'); % 绘制PTS处理后的PAPR CCDF
legend('Orignal','PTS'); % 图例
title('V=4'); % 标题
xlabel('PAPR0 [dB]'); % X轴标签
ylabel('CCDF (Pr[PAPR>PAPR0])'); % Y轴标签
grid on; % 显示网格



if Lfft==128
   save QP_PTS1.mat PAPR1 PAPR2 cdf1 cdf2
end
if Lfft==256
   save QP_PTS2.mat PAPR1 PAPR2 cdf1 cdf2
end
if Lfft==512
   save QP_PTS3.mat PAPR1 PAPR2 cdf1 cdf2
end





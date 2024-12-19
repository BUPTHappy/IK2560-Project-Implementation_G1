clc;                      % ��������
clear;                    % ������б���
close all;                % �ر�����ͼ�δ���
warning off;              % �رվ�����Ϣ
addpath(genpath(pwd));    % ����ǰĿ¼��ӵ�����·��
rng('default');           % ����������������ĳ�ʼ״̬ΪĬ��ֵ
rng(1);  

% FFT��С����
Lfft_values = [64, 128, 256, 512, 1024];
V = 4;            % ѡ�������
MAP_qpsk  = [1 -1 j -j];  % QPSK���Ƽ�
Phases    = [1 -1];       % ��λ���� W = 2

Nbits     = 10000;      % ����������

% ���������и��� = W ^(V-1) = 8
Pchos     = [1 1 1 1; 1 1 1 2; 1 1 2 1; 1 2 1 1; 2 1 1 1;...
             1 1 2 2; 1 2 1 2; 1 2 2 1; 2 2 1 1; 2 1 2 1; 2 1 1 2;...
             2 2 2 1; 2 2 1 2; 2 1 2 2; 1 2 2 2; 2 2 2 2];
Lchos     = 16;         % ѡ�����еĳ���

% ����ͼ��
figure;
hold on;

% Ԥ������ɫ
colors = lines(length(Lfft_values));  % ʹ��MATLAB��'lines'��ɫ��

% �洢ͼ���ַ���
legend_entries = {};

% ѭ��������ͬ��FFT��С
for idx = 1:length(Lfft_values)
    Lfft = Lfft_values(idx);  % ��ǰFFT��С
    papr0     = zeros(1,Nbits); % �洢ԭʼPAPRֵ
    papr_pts  = zeros(1,Nbits); % �洢PTS������PAPRֵ

    for jj=1:Nbits
        Index        = floor(length(MAP_qpsk)*rand(1,Lfft))+1; % ���ѡ��QPSK��������
        X            = MAP_qpsk(Index(1,:)); % ԭʼƵ���ź�
        x            = ifft(X,[],2); % ʱ���ź�
        Pow1         = abs(x.^2); % �źŹ���
        Pow2         = max(Pow1,[],2); % ��ֵ����
        Pow3         = mean(Pow1,2);   % ƽ������
        papr0(jj)    = 10*log10(Pow2./Pow3); % ����ԭʼPAPR

        % PTS�ָ�
        A = zeros(V,Lfft); % ��ʼ��A����
        for v=1:V
            A(v,(1+(v-1)*Lfft/4):1:v*Lfft/4) = X((1+(v-1)*Lfft/4):1:v*Lfft/4); % ��X�ֳ�V�����ַ���A��
        end
        a    = ifft(A,[],2); % ��ÿ��������IFFT
        % ���޺���
        Vmin = 10; % ��ʼ����Сֵ
        % Ѱ�����Ÿ�����Ϣ
        for n=1:Lchos
            temp_phase = Phases(Pchos(n,:)).'; % ��ȡ��λ���
            b = repmat(temp_phase,1,Lfft); % ������λ���
            c = a.*b; % ���
            d = sum(c); % ���
            e = abs(d); % ȡ����ֵ
            temp_max = max(e); % �������ֵ
            if temp_max<Vmin
                Vmin = temp_max; % ������Сֵ
                Best_n = n; % ��¼���ѡ��
            end
        end
        aa = sum(a.*repmat(Phases(Pchos(Best_n,:)).',1,Lfft)); % ���ѡ����λ���������

        Pow1 = abs(aa.^2); % �µ��źŹ���
        Pow2   = max(Pow1,[],2); % �µķ�ֵ����
        Pow3   = mean(Pow1,2);   % �µ�ƽ������
        papr_pts(jj) = 10*log10(Pow2./Pow3); % �����µ�PAPR
    end

    % �����ۻ��ֲ�����
    [cdf1, PAPR1] = ecdf(papr0); % ԭʼPAPR��CDF
    [cdf2, PAPR2] = ecdf(papr_pts);     % PTS������PAPR CDF

    % ���Ʋ�ͬFFT��С�Ľ��
    semilogy(PAPR1(1:end), 1-cdf1(1:end), 'Color', colors(idx,:), 'LineWidth', 2);
    semilogy(PAPR2(1:end), 1-cdf2(1:end), '--', 'Color', colors(idx,:), 'LineWidth', 2);  

    % ���ͼ����Ŀ
    legend_entries{end+1} = sprintf('FFT=%d (Original)', Lfft);
    legend_entries{end+1} = sprintf('FFT=%d (PTS)', Lfft);
    
end

% ����ͼ��
legend(legend_entries, 'Location', 'Best');
title('Comparison of PAPR for different FFT sizes with and without PTS');
xlabel('PAPR [dB]');
ylabel('CCDF (Pr[PAPR>PAPR0])');
grid on;
hold off;

% ���Ʋ�ͬFFT��С��PAPRֱ��ͼ
figure;
hold on;

% ����͸����
alpha_value = 0.5; % ͸��������Ϊ50%

% ��ɫ����ɫ����ɫ
color_original = [1, 1, 0];  % ��ɫ
color_pts = [0, 0, 1];       % ��ɫ

% ʹ�ò�ͬ��ɫ����ԭʼPAPR��PTS������PAPR
for idx = 1:length(Lfft_values)
    Lfft = Lfft_values(idx);
    
    % ԭʼPAPR��ֱ��ͼ����ɫ��
    histogram(papr0, 'Normalization', 'pdf', 'EdgeColor', 'none', ...
        'FaceColor', color_original, 'BinWidth', 0.1, 'FaceAlpha', alpha_value);
    
    % PTS������PAPR��ֱ��ͼ����ɫ��
    histogram(papr_pts, 'Normalization', 'pdf', 'EdgeColor', 'none', ...
        'FaceColor', color_pts, 'BinWidth', 0.1, 'LineStyle', '--', 'FaceAlpha', alpha_value);
end

% ����ֶ����õ�ͼ��
legend_entries = {'Original PAPR (Yellow)', 'PTS PAPR (Blue)'};

% ����ͼ��
legend(legend_entries, 'Location', 'Best');
title('PAPR Histogram for Different FFT Sizes with and without PTS');
xlabel('PAPR [dB]');
ylabel('Probability Density');
grid on;
hold off;


% ����ͼ�δ���
figure;
hold on;
% ʹ��subplot���ƶ����ͼ
for idx = 1:length(Lfft_values)
    Lfft = Lfft_values(idx);  % ��ǰFFT��С
    
    % ���ѡ��QPSK��������
    Index = floor(length(MAP_qpsk)*rand(1, Lfft)) + 1;
    X = MAP_qpsk(Index(1,:)); % Ƶ���ź�
    
    % ͨ��IFFT����ʱ���ź�
    x = ifft(X, Lfft);
    
    % ������ͼ
    subplot(2, 3, idx);  % 2��3�е���ͼ����idx����ͼ
    plot(real(x), 'b-', 'LineWidth', 2);
    hold on;
    plot(imag(x), 'r--', 'LineWidth', 2);
    title(['OFDM Signal Time-domain (FFT=' num2str(Lfft) ')']);
    xlabel('Sample Index');
    ylabel('Amplitude');
    legend({'Real Part', 'Imaginary Part'}, 'Location', 'Best');
    grid on;
end

% ��������ͼ��
sgtitle('OFDM Signal Time-domain for Different FFT Sizes');
hold off;
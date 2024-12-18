clc;                      % ��������
clear;                    % ������б���
close all;                % �ر�����ͼ�δ���
warning off;              % �رվ�����Ϣ
addpath(genpath(pwd));    % ����ǰĿ¼��ӵ�����·��
rng('default');           % ����������������ĳ�ʼ״̬ΪĬ��ֵ
rng(1);  


Lfft      = 128;          % FFT�Ĵ�С
V         = 4;            % ѡ�������
MAP_qpsk  = [1 -1 j -j];  % QPSK���Ƽ�
Phases    = [1 -1];       % ��λ���� W = 2

Nbits     = 10000;      % ����������

% ���������и��� = W ^(V-1) = 8
% ���QPSK����16����λ�仯���
Pchos     = [1 1 1 1; 1 1 1 2; 1 1 2 1; 1 2 1 1; 2 1 1 1;...
             1 1 2 2; 1 2 1 2; 1 2 2 1; 2 2 1 1; 2 1 2 1; 2 1 1 2;...
             2 2 2 1; 2 2 1 2; 2 1 2 2; 1 2 2 2; 2 2 2 2];
Lchos     = 16;         % ѡ�����еĳ���


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
    % ���޺������˴�ʡ�ԣ�
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

% ��ͼ
semilogy(PAPR1(1:1:end),1-cdf1(1:1:end),'-b',... % ����ԭʼPAPR��CCDF
         PAPR2(1:1:end),1-cdf2(1:1:end),'-r'); % ����PTS������PAPR CCDF
legend('Orignal','PTS'); % ͼ��
title('V=4'); % ����
xlabel('PAPR0 [dB]'); % X���ǩ
ylabel('CCDF (Pr[PAPR>PAPR0])'); % Y���ǩ
grid on; % ��ʾ����



if Lfft==128
   save QP_PTS1.mat PAPR1 PAPR2 cdf1 cdf2
end
if Lfft==256
   save QP_PTS2.mat PAPR1 PAPR2 cdf1 cdf2
end
if Lfft==512
   save QP_PTS3.mat PAPR1 PAPR2 cdf1 cdf2
end





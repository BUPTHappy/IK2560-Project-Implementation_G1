clear; clc;
%% Parameters
Nffts=2.^[6:10];
Npsk=2;
M=2^Npsk;
Nblk=1e3;
Z_dBs=4:0.1:10;
N_Z_dBs=length(Z_dBs);
CCDF_formula=@(N,sigma2,z) 1-((1-exp(-z.^2/(2*sigma2))).^N);

%% main function
for kk=1:length(Nffts)
    Nfft=Nffts(kk);
    x=zeros(Nfft,Nblk);
    x_CF=zeros(1,Nblk);
    CCDF_simulation=zeros(1,N_Z_dBs);
    for ii=1:Nblk
        X_mod=ModSymbolGenrator(Npsk,Nfft);
        x(:,ii)=ifft(ifftshift(X_mod),Nfft) .*sqrt(Nfft);
        x_CF(ii)=PAPR_dB(x(:,ii));
    end
    sigma2=mean(mean(abs(x)))^2/(pi/2);
    CCDF_theoretical=CCDF_formula(Nfft,sigma2,10.^(Z_dBs/20));
    for jj= 1:N_Z_dBs
        CCDF_simulation(jj)=sum(x_CF>Z_dBs(jj))/Nblk;
    end
    semilogy(Z_dBs,CCDF_theoretical,'-');
    hold on;
    grid on;
    semilogy(Z_dBs,CCDF_simulation,'*');
end
title('OFDM system with N-point FFT');
xlabel('PAPR0[dB]');
ylabel('CCDF=Probability(PAPR>PAPR0)');
legend('Theoretical','Simulation');
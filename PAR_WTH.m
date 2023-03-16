clc;
clear;

load('OFDM_IM_Para.mat');

%% 迭代参数
iter = 1;  % # Iterations 迭代次数
PAR0_dB = 0:30;
len =2.5;
OFDM_IM_Para.N_sc_used = 22;
CCDF_k =  zeros(4,size(PAR0_dB,2));
CCDF_WTH_k =  zeros(4,size(PAR0_dB,2));
for k =20:2:20
    OFDM_IM_Para.M_QAM = 6;
    OFDM_IM_Para.K_sc_used = k;
    OFDM_IM_Para.Bit_Loading_SC = zeros(1, OFDM_IM_Para.N_fft);
    OFDM_IM_Para.Bit_Loading_SC_active = OFDM_IM_Para.Bit_Loading_SC;
    OFDM_IM_Para.Bit_Loading_SC(1: OFDM_IM_Para.N_sc_used) = OFDM_IM_Para.M_QAM * ones(1, OFDM_IM_Para.N_sc_used);
    OFDM_IM_Para.Bit_Loading_SC = circshift(OFDM_IM_Para.Bit_Loading_SC, [0, floor(-OFDM_IM_Para.N_sc_used/2)]);
    OFDM_IM_Para.Bit_Loading_SC_active(1: OFDM_IM_Para.K_sc_used) = OFDM_IM_Para.M_QAM * ones(1, OFDM_IM_Para.K_sc_used);
    OFDM_IM_Para.Bit_Loading_SC_active= circshift(OFDM_IM_Para.Bit_Loading_SC_active, [0, floor(-OFDM_IM_Para.K_sc_used/2)]);
    OFDM_IM_Para.Bit_index = sort((1:OFDM_IM_Para.N_fft) - (OFDM_IM_Para.N_fft - OFDM_IM_Para.N_sc_used),"descend");
    OFDM_IM_Para.Bit_index = circshift(OFDM_IM_Para.Bit_index, [0, floor(-OFDM_IM_Para.N_sc_used/2)]);

    OFDM_IM_Para.N_Bit_Symbol_info = sum(OFDM_IM_Para.Bit_Loading_SC_active);         %number of transmitted information bits in each OFDM IM symbol
    OFDM_IM_Para.N_Bit_Symbol_index = floor(log2(nchoosek(OFDM_IM_Para.N_sc_used,OFDM_IM_Para.K_sc_used)));        %number of transmitted index bits in each OFDM IM symbol
    OFDM_IM_Para.N_Bit_Symbol = OFDM_IM_Para.N_Bit_Symbol_info + OFDM_IM_Para.N_Bit_Symbol_index;                % number of transmitted bits in each OFDM IM symbol每个OFDM符号传输的比特数量

    CCDF = zeros(1,size(PAR0_dB,2));
    CCDF_WTH = zeros(1,size(PAR0_dB,2));
    Double_CCDF_WTH = zeros(1,size(PAR0_dB,2));

    for s1 = PAR0_dB + 1
        ccdf_iter = zeros(1,iter);
        ccdf_wth_iter = zeros(1,iter);
        for s2 = 1:iter
            fprintf('== PAR0(dB) is %g and iteration is %g == \n',PAR0_dB(s1)./len,s2)
            [Tx_Signal_OS,Tx_Signal1,Tx_Signal2] = PAR_OFDM_IM_TX(OFDM_IM_Para,9);
            PARdB = PAR_caculate(Tx_Signal_OS,OFDM_IM_Para.N_fft);
            PARdB_WTH = PAR_caculate(Tx_Signal1,OFDM_IM_Para.N_fft);
            PARdB_Double_WTH = PAR_caculate(Tx_Signal2,OFDM_IM_Para.N_fft);
            ccdf_iter(s2) = sum(sum(PARdB > PAR0_dB(s1)./len))/(OFDM_IM_Para.N_Symbol);
            ccdf_wth_iter(s2) = sum(sum(PARdB_WTH> PAR0_dB(s1)./len))/(OFDM_IM_Para.N_Symbol);
            double_ccdf_wth_iter(s2) = sum(sum(PARdB_Double_WTH> PAR0_dB(s1)./len))/(OFDM_IM_Para.N_Symbol);
        end
        CCDF(s1) = sum(ccdf_iter)/iter;
        CCDF_WTH(s1) = sum(ccdf_wth_iter)/iter;
        Double_CCDF_WTH(s1) = sum(double_ccdf_wth_iter)/iter;       
    end
CCDF_k((k-12)/2,:) = CCDF;
CCDF_WTH_k((k-12)/2,:) = CCDF_WTH;
Double_CCDF_WTH_k((k-12)/2,:) = Double_CCDF_WTH;
end    
    
    % [cdf0, PAPR0] = ecdf(PARdB);%计算CCDF
    % [cdf1, PAPR1] = ecdf(PARdB_WTH);%计算CCDF


figure(1);
% semilogy(PAPR0,1-cdf0,'r*-','LineWidth',1.5,'MarkerSize',12)
semilogy(PAR0_dB./len,CCDF_k(1,:),'r.-','LineWidth',1.5,'MarkerSize',12)
hold on
% semilogy(PAPR1,1-cdf1,'ks-','LineWidth',1.5,'MarkerSize',12) 
semilogy(PAR0_dB./len,CCDF_WTH_k(1,:),'k.-','LineWidth',1.5,'MarkerSize',12)
hold on

semilogy(PAR0_dB./len,CCDF_k(2,:),'b.-','LineWidth',1.5,'MarkerSize',12)
hold on
% semilogy(PAPR1,1-cdf1,'ks-','LineWidth',1.5,'MarkerSize',12) 
semilogy(PAR0_dB./len,CCDF_WTH_k(2,:),'g.-','LineWidth',1.5,'MarkerSize',12)
hold on

semilogy(PAR0_dB./len,CCDF_k(3,:),'m.-','LineWidth',1.5,'MarkerSize',12)
hold on
% semilogy(PAPR1,1-cdf1,'ks-','LineWidth',1.5,'MarkerSize',12) 
semilogy(PAR0_dB./len,CCDF_WTH_k(3,:),'y.-','LineWidth',1.5,'MarkerSize',12)
hold on

semilogy(PAR0_dB./len,CCDF_k(4,:),'.-','LineWidth',1.5,'MarkerSize',12)
hold on
% semilogy(PAPR1,1-cdf1,'ks-','LineWidth',1.5,'MarkerSize',12) 
semilogy(PAR0_dB./len,CCDF_WTH_k(4,:),'.-','LineWidth',1.5,'MarkerSize',12)
hold on

grid on
title("OFDM-IM")
xlabel('PAR(dB)')
ylabel('CCDF')
legend('N=22,K=14','WTH N=22,K=14','N=22,K=16','WTH N=22,K=16','N=22,K=18','WTH N=22,K=18','N=22,K=20','WTH N=22,K=20')
legend('N=22,K=14','Double WTH N=22,K=14','N=22,K=16','Double WTH N=22,K=16','N=22,K=18','Double WTH N=22,K=18','N=22,K=20','Double WTH N=22,K=20')
legend('WTH N=22,K=14','Double WTH N=22,K=14','WTH N=22,K=16','Double WTH N=22,K=16','WTH N=22,K=18','Double WTH N=22,K=18','WTH N=22,K=20','Double WTH N=22,K=20')
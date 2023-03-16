clc;
clear;

load('OFDM_IM_Para.mat');
[Tx_Signal_OS,Tx_Signal1,Tx_Signal2] = Double_WHT_OFDM_IM_TX(OFDM_IM_Para,9);
load('index_bit.mat');
load('index_sym.mat');
load('info_bit.mat');
load('tx_symbol.mat');


%% 迭代参数
iter = 1;  % # Iterations 迭代次数
EbN0dB = 0:1:25; %信噪比（dB）
iter = 1;  % # Iterations 迭代次数
% EbN0dB = 10:5:15; %信噪比（dB）
EsN0 = 10.^(EbN0dB/10); %信噪比
sigma = sqrt(sum(abs(Tx_Signal_OS).^2)./EsN0./size(Tx_Signal_OS,2)); % additive noise variance %加性白噪声的方差

%% ==================== Loop for SNR =========================
PEP = zeros(1,size(sigma,2)); % index symbol error IEP 索引符号误差 IEP
OFDM_SER = zeros(1,size(sigma,2)); % M-ary symbol error  OFDM的误码率
Total_SER = zeros(1,size(sigma,2)); % SEP overall 总误码率
BER=zeros(1,size(sigma,2));
BER1=zeros(1,size(sigma,2)); % index bit error rate 索引误比特率
BER2=zeros(1,size(sigma,2)); % M-ary bit error rate 调制误比特率
PEP_ML = zeros(1,size(sigma,2)); % index symbol error IEP 索引符号误差 IEP
OFDM_SER_ML = zeros(1,size(sigma,2)); % M-ary symbol error  OFDM的误码率
Total_SER_ML = zeros(1,size(sigma,2)); % SEP overall 总误码率
BER_ML=zeros(1,size(sigma,2));
BER1_ML=zeros(1,size(sigma,2)); % index bit error rate 索引误比特率
BER2_ML=zeros(1,size(sigma,2)); % M-ary bit error rate 调制误比特率
PEP_LLR = zeros(1,size(sigma,2)); % index symbol error IEP 索引符号误差 IEP
OFDM_SER_LLR = zeros(1,size(sigma,2)); % M-ary symbol error  OFDM的误码率
Total_SER_LLR = zeros(1,size(sigma,2)); % SEP overall 总误码率
BER_LLR=zeros(1,size(sigma,2));
BER1_LLR=zeros(1,size(sigma,2)); % index bit error rate 索引误比特率
BER2_LLR=zeros(1,size(sigma,2)); % M-ary bit error rate 调制误比特率

t1 = tic;
for s1 = 1:size(sigma,2)
    fprintf('== EbN0(dB) is %g == \n',EbN0dB(s1))
        %% ==================== Loop for iteration 初始化 =======================
    symerr_im = zeros(1,iter);
    symerr_ofdm = zeros(1,iter);    %iter迭代次数
    symerr_iter= zeros(1,iter);
    BER_iter= zeros(1,iter);
    BER_iter_1= zeros(1,iter);
    BER_iter_2= zeros(1,iter);
    symerr_im_ML = zeros(1,iter);
    symerr_ofdm_ML = zeros(1,iter);    %iter迭代次数
    symerr_iter_ML= zeros(1,iter);
    BER_iter_ML= zeros(1,iter);
    BER_iter_1_ML= zeros(1,iter);
    BER_iter_2_ML= zeros(1,iter);
    symerr_im_LLR = zeros(1,iter);
    symerr_ofdm_LLR = zeros(1,iter);    %iter迭代次数
    symerr_iter_LLR= zeros(1,iter);
    BER_iter_LLR= zeros(1,iter);
    BER_iter_1_LLR= zeros(1,iter);
    BER_iter_2_LLR= zeros(1,iter);
    for s2 = 1:iter
        fprintf('== EbN0(dB) is %g and iteration is %g == \n',EbN0dB(s1),s2)
        noise = sigma(s1)*(randn(size(Tx_Signal_OS))+1i*randn(size(Tx_Signal_OS)));
        y = Tx_Signal_OS + noise;
        y1 = Tx_Signal1 + noise;
        y2 = Tx_Signal2 + noise;
        OFDM_IM_Para.detection = 3;
        [Rx_bin_index ,Rx_bin_info,index_symbol_de,infor_symbol_de] = OFDM_IM_RX(OFDM_IM_Para, y,sigma(s1)^2,9,EbN0dB(s1));
        [Rx_bin_index_ML ,Rx_bin_info_ML,index_symbol_de_ML,infor_symbol_de_ML] = WHT_OFDM_IM_RX(OFDM_IM_Para, y1,sigma(s1)^2,9,EbN0dB(s1));
        [Rx_bin_index_LLR ,Rx_bin_info_LLR,index_symbol_de_LLR,infor_symbol_de_LLR] = Double_WHT_OFDM_IM_RX(OFDM_IM_Para, y2,sigma(s1)^2,9,EbN0dB(s1));
%         OFDM_IM_Para.detection = 2;
%         [Rx_bin_index_LLR ,Rx_bin_info_LLR,index_symbol_de_LLR,infor_symbol_de_LLR] = OFDM_IM_RX(OFDM_IM_Para, y,sigma(s1)^2,9,EbN0dB(s1));
%         
%         OFDM_IM_Para.detection = 1;
%         [Rx_bin_index_ML ,Rx_bin_info_ML,index_symbol_de_ML,infor_symbol_de_ML] = OFDM_IM_RX(OFDM_IM_Para, y,sigma(s1)^2,9,EbN0dB(s1));
%         
        % QAM调制符号差错个数
        ofdm_symerr = sum(sum(tx_symbol~=infor_symbol_de));
        % 索引符号差错个数
        ind_symerr = sum(index_sym~=index_symbol_de);
        %QAM调制比特差错个数
        info_bit_err=sum(sum(info_bit~=Rx_bin_info));
        %索引比特差错个数
        index_bit_err=sum(sum(index_bit~=Rx_bin_index)); 
        
        %SER
        symerr_im(s2) = ind_symerr/OFDM_IM_Para.N_Symbol;
        symerr_ofdm(s2) = ofdm_symerr/(OFDM_IM_Para.K_sc_used*OFDM_IM_Para.N_Symbol);
        symerr_iter(s2) = (ind_symerr+ofdm_symerr)/(OFDM_IM_Para.N_Symbol+OFDM_IM_Para.K_sc_used*OFDM_IM_Para.N_Symbol);
        
        %BER
        BER_iter(s2)=(info_bit_err+index_bit_err)./(OFDM_IM_Para.N_Bit_Symbol*OFDM_IM_Para.N_Symbol);
        BER_iter_1(s2) = index_bit_err./(OFDM_IM_Para.N_Bit_Symbol_index*OFDM_IM_Para.N_Symbol);
        BER_iter_2(s2) = info_bit_err./(OFDM_IM_Para.N_Bit_Symbol_info*OFDM_IM_Para.N_Symbol);
    
        % QAM调制符号差错个数
        ofdm_symerr_ML = sum(sum(tx_symbol~=infor_symbol_de_ML));
        % 索引符号差错个数
        ind_symerr_ML = sum(index_sym~=index_symbol_de_ML);
        %QAM调制比特差错个数
        info_bit_err_ML=sum(sum(info_bit~=Rx_bin_info_ML));
        %索引比特差错个数
        index_bit_err_ML=sum(sum(index_bit~=Rx_bin_index_ML)); 
        
        %SER
        symerr_im_ML(s2) = ind_symerr_ML/OFDM_IM_Para.N_Symbol;
        symerr_ofdm_ML(s2) = ofdm_symerr_ML/(OFDM_IM_Para.K_sc_used*OFDM_IM_Para.N_Symbol);
        symerr_iter_ML(s2) = (ind_symerr_ML+ofdm_symerr_ML)/(OFDM_IM_Para.N_Symbol+OFDM_IM_Para.K_sc_used*OFDM_IM_Para.N_Symbol);
        
        %BER
        BER_iter_ML(s2)=(info_bit_err_ML+index_bit_err_ML)./(OFDM_IM_Para.N_Bit_Symbol*OFDM_IM_Para.N_Symbol);
        BER_iter_1_ML(s2) = index_bit_err_ML./(OFDM_IM_Para.N_Bit_Symbol_index*OFDM_IM_Para.N_Symbol);
        BER_iter_2_ML(s2) = info_bit_err_ML./(OFDM_IM_Para.N_Bit_Symbol_info*OFDM_IM_Para.N_Symbol);        

        % QAM调制符号差错个数
        ofdm_symerr_LLR = sum(sum(tx_symbol~=infor_symbol_de_LLR));
        % 索引符号差错个数
        ind_symerr_LLR = sum(index_sym~=index_symbol_de_LLR);
        %QAM调制比特差错个数
        info_bit_err_LLR=sum(sum(info_bit~=Rx_bin_info_LLR));
        %索引比特差错个数
        index_bit_err_LLR=sum(sum(index_bit~=Rx_bin_index_LLR)); 
        
        %SER
        symerr_im_LLR(s2) = ind_symerr_LLR/OFDM_IM_Para.N_Symbol;
        symerr_ofdm_LLR(s2) = ofdm_symerr_LLR/(OFDM_IM_Para.K_sc_used*OFDM_IM_Para.N_Symbol);
        symerr_iter_LLR(s2) = (ind_symerr_LLR+ofdm_symerr_LLR)/(OFDM_IM_Para.N_Symbol+OFDM_IM_Para.K_sc_used*OFDM_IM_Para.N_Symbol);
        
        %BER
        BER_iter_LLR(s2)=(info_bit_err_LLR+index_bit_err_LLR)./(OFDM_IM_Para.N_Bit_Symbol*OFDM_IM_Para.N_Symbol);
        BER_iter_1_LLR(s2) = index_bit_err_LLR./(OFDM_IM_Para.N_Bit_Symbol_index*OFDM_IM_Para.N_Symbol);
        BER_iter_2_LLR(s2) = info_bit_err_LLR./(OFDM_IM_Para.N_Bit_Symbol_info*OFDM_IM_Para.N_Symbol);
    
    end
    %% =============average bit error rate================
    PEP(s1) = sum(symerr_im)/iter;
    OFDM_SER(s1) = sum(symerr_ofdm)/iter;
    Total_SER(s1) = sum(symerr_iter)/iter;
    BER(s1)= sum(BER_iter)./iter;
    BER1(s1)= sum(BER_iter_1)./iter;
    BER2(s1)= sum(BER_iter_2)./iter;
    PEP_ML(s1) = sum(symerr_im_ML)/iter;
    OFDM_SER_ML(s1) = sum(symerr_ofdm_ML)/iter;
    Total_SER_ML(s1) = sum(symerr_iter_ML)/iter;
    BER_ML(s1)= sum(BER_iter_ML)./iter;
    BER1_ML(s1)= sum(BER_iter_1_ML)./iter;
    BER2_ML(s1)= sum(BER_iter_2_ML)./iter;
    PEP_LLR(s1) = sum(symerr_im_ML)/iter;
    OFDM_SER_LLR(s1) = sum(symerr_ofdm_LLR)/iter;
    Total_SER_LLR(s1) = sum(symerr_iter_LLR)/iter;
    BER_LLR(s1)= sum(BER_iter_LLR)./iter;
    BER1_LLR(s1)= sum(BER_iter_1_LLR)./iter;
    BER2_LLR(s1)= sum(BER_iter_2_LLR)./iter;
end
 toc(t1)
figure;
semilogy(EbN0dB,BER,' rs-','LineWidth',1.5,'MarkerSize',12)
hold on
semilogy(EbN0dB,Total_SER,' b+-','LineWidth',1.5,'MarkerSize',12)
% hold on
% semilogy(EbN0dB,BER2,'m d-','LineWidth',1.5,'MarkerSize',12)
% hold on
% semilogy(EbN0dB,BER,'k s-','LineWidth',1.5,'MarkerSize',12)
% hold on
% semilogy(EbN0dB,Total_SER,'m o-','LineWidth',1.5,'MarkerSize',12)
% hold on
semilogy(EbN0dB,BER_ML,' ks-','LineWidth',1.5,'MarkerSize',12)
hold on
semilogy(EbN0dB,Total_SER_ML,' mo-','LineWidth',1.5,'MarkerSize',12)
hold on
% semilogy(EbN0dB,BER2_ML,'c h-','LineWidth',1.5,'MarkerSize',12)
% hold on


semilogy(EbN0dB,BER_LLR,'g d-','LineWidth',1.5,'MarkerSize',12)
hold on
semilogy(EbN0dB,Total_SER_LLR,'c h-','LineWidth',1.5,'MarkerSize',12)
hold on
% axis([0 27 10^-8 10^0])
grid on
hold on
title('')
xlabel('Es/No (dB)')
ylabel('BER/SER')
% legend('MSP BER','MSP SER','MSP no power threshold  BER','MSP no power threshold SER')
% legend('MSP no power threshold  BER','MSP no power threshold SER','MSP BER','MSP SER')
% legend('MSP BER','MSP SER','ML BER','ML SER','LLR BER','LR')
legend('NO WHT BER','NO WHT SER','WHT BER','WHT SER','Double WHT BER','Double WHT SER')
% legend('N=16,K=8','WHT N=16,K=8','N=16,K=10','WHT NLR SE=16,K=10','N=16,K=12','WHT N=16,K=12','N=16,K=14','WHT N=16,K=14')
% title("OFDM-IM")
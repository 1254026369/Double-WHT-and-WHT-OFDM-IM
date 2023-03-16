function [Rx_bin_index, Rx_bin_info,index_symbol_de,infor_symbol_de] = Double_WHT_OFDM_IM_RX(OFDM_IM_Para, Rx_Sample,N_sigma,bit,EbN0dB)
%% OFDM_IM_RX  OFDM IM接收机

%% ---------- OFDM Parameters
N_Symbol = OFDM_IM_Para.N_Symbol;
detection = OFDM_IM_Para.detection;
N_frm = OFDM_IM_Para.N_frm;
Nb_fft = OFDM_IM_Para.Nb_fft;
CP = OFDM_IM_Para.CP;
CPshift = OFDM_IM_Para.CPshift;                          % CP offset for OFDM demodulation
ClippingRatio = OFDM_IM_Para.ClippingRatio;
N_sc_used = OFDM_IM_Para.N_sc_used;
K_sc_used = OFDM_IM_Para.K_sc_used;
N_fft = OFDM_IM_Para.N_fft;
Bit_Loading_SC = OFDM_IM_Para.Bit_Loading_SC;
Bit_Loading_SC_active = OFDM_IM_Para.Bit_Loading_SC_active;
N_Bit_Symbol_info = OFDM_IM_Para.N_Bit_Symbol_info;
N_Bit_Symbol_index = OFDM_IM_Para.N_Bit_Symbol_index;
N_Bit_Symbol = OFDM_IM_Para.N_Bit_Symbol;
BW_OFDM = OFDM_IM_Para.BW_OFDM;
Clipping_Agg = OFDM_IM_Para.Clipping_Agg;
Bit_index = OFDM_IM_Para.Bit_index;
Quantizationbit = OFDM_IM_Para.Quantizationbit;

%% 组合法产生索引映射
index_all = Combin_Md(N_sc_used,K_sc_used);     % Combin_Md函数索引映射方式的实现

%% 产生检测符号
for qq=1:2^Bit_Loading_SC(1)
    sym_test(qq)=qammod(qq-1,2^Bit_Loading_SC(1),'gray');
end

% % 功率归一化
% for k= 1: 2^Bit_Loading_SC(1)
%     Mod_Bit = Bit_Loading_SC(k);
%     switch Mod_Bit
%         case 0
%             P_mean = 1;
%         case 2
%             P_mean = 2;
%         case 4
%             P_mean = 10;     % normalized powers for different QAM modulation formats归一化功率
%         case 5
%             P_mean = 20;
%         case 6
%             P_mean = 42;
%         case 7
%             P_mean = 82;
%         case 8
%             P_mean = 170;
%         otherwise
%             error('Incorrect Subcarrier Modulation Format!')
%     end
%     sym_test = sym_test/sqrt(P_mean);
% end

%% --------------------ADC--------------------------%
N_bits = Quantizationbit; % 8-bits DAC/ADC
ss = Rx_Sample;
ss = ss - mean(ss);
ss_i = imag(ss);
ss_r = real(ss);
V_max_r= max(abs(ss_r));
V_max_i = max(abs(ss_i));
% V_max = max(abs(ss));
% d_V = 2 * V_max / 2^N_bits;
d_V_r = 2 * V_max_r / 2^N_bits;
d_V_i= 2 * V_max_i / 2^N_bits;
Rx_Sample_r = (floor(ss_r/d_V_r) + 1/2) * d_V_r; % ADC
Rx_Sample_i = (floor(ss_i/d_V_i) + 1/2) * d_V_i; % ADC
Rx_Sample = Rx_Sample_r+Rx_Sample_i.*j;


%% 去除循环前缀 + FFT
 Rx_Sym_FFT = zeros(N_Symbol, N_fft, N_frm);
 NN = 2.^ceil(log2(N_sc_used));
Hada = hadamard(NN);
 
for k = 1: N_frm
    Rx_Temp = vec2mat(Rx_Sample(k, :), N_fft + ceil(N_fft * CP));
    Rx_Temp = Rx_Temp(:,ceil(N_fft * CP) + 1 - CPshift: N_fft + ceil(N_fft * CP) - CPshift); % cut the CP

    Rx_Temp =  Rx_Temp * Hada; 
    OFDM_Rx_LTE = fft(Rx_Temp, N_fft, 2)/N_fft;
    Rx_Sym_FFT(:, :, k) = OFDM_Rx_LTE;    
end
% 
% Hada = hadamard(N_sc_used);
Rx_Sym_FFT =  Rx_Sym_FFT * Hada; 
% if NN~= N_sc_used
%     Rx_Sym_FFT = rx_Sym_FFT(:,1:N_sc_used);
% end


%% 信道估计
%----------------- Channel Estimation-------------------%
load('Tx_Sym_FFT.mat');

Rx_Sym_SC_Eq = zeros(size(Rx_Sym_FFT, 1), N_sc_used, N_frm);
Rx_figure = zeros(size(Rx_Sym_FFT, 1), N_sc_used, N_frm);

for k= 1: N_frm    
    Rx_FFT = Rx_Sym_FFT(:, :, k);
    Rx_FFT = circshift(Rx_FFT, [0, ceil(N_sc_used/2)]);
    Rx = Rx_FFT(:, 1: N_sc_used);
    %figure(); plot(Rx(:), '.');
    Rx_figure(:, :, k) =  Rx;
    Tx_FFT = Tx_Sym_FFT(:, :, k);
    Tx_FFT = circshift(Tx_FFT, [0, ceil(N_sc_used/2)]);
    Tx = Tx_FFT(:, 1: N_sc_used);  
    mm = Rx./Tx;
    mm(find(isnan(mm)==1)) = 0;
    mm(find(isinf(mm)==1))=0;
    CRdata = sum(mm,1)./sum((mm~=0),1);
    Rx_Eq = Rx.*repmat(conj(CRdata)./abs(CRdata).^2, size(Rx, 1), 1);
    Rx_Sym_SC_Eq(:, :, k) = Rx_Eq;   
end
% 
% scatterplot(Rx_Sym_SC_Eq(:));
% scatterplot(Tx(:));

%% 接收机检测IM检测（两种判决方式）
h = 1;
Mod_Bit = Bit_Loading_SC(1);
%MLtest
% for k = 1: N_frm
%     [index_symbol_de(:,:,k),re_sym] = ML_detection(K_sc_used,N_Bit_Symbol_index,index_all,Tx,h,2^Bit_Loading_SC(1),N_Symbol,sym_test);
%     infor_symbol_de(:,:,k) = qamdemod(re_sym, 2^Mod_Bit, 'gray');%格雷编码
% end

% %LLR test
% N_sigma = 1;
% for k = 1: N_frm
%     [index_symbol_de(:,:,k),infor_symbol_de(:,:,k)] = LLR_detection(K_sc_used,Tx,h,2^Bit_Loading_SC(1),N_Symbol,sym_test,N_sigma);
% end


 if detection == 1      %ML判决
    for k = 1: N_frm
        [index_symbol_de(:,:,k),re_sym] = ML_detection(K_sc_used,N_Bit_Symbol_index,index_all,Rx_Sym_SC_Eq(:, :, k),h,2^Bit_Loading_SC(1),N_Symbol,sym_test);
        infor_symbol_de(:,:,k) = qamdemod(re_sym, 2^Mod_Bit, 'gray');%格雷编码
    end
 elseif detection == 2         %LLR判决
%     N_sigma = 1; 测试
    for k = 1: N_frm
        [index_symbol_de(:,:,k),infor_symbol_de(:,:,k)] = LLR_detection(K_sc_used,Rx_Sym_SC_Eq(:, :, k),h,2^Bit_Loading_SC(1),N_Symbol,sym_test,abs(mean(conj(CRdata)./abs(CRdata).^2))^2*N_sigma);
    end
 elseif detection == 3         %MP判决
    for k = 1: N_frm
        [index_symbol_de(:,:,k),infor_symbol_de(:,:,k)] = MP_detection(K_sc_used,Rx_Sym_SC_Eq(:, :, k),2^Bit_Loading_SC(1),N_Symbol,abs(mean(conj(CRdata)./abs(CRdata).^2))^2*N_sigma,EbN0dB);
        %[index_symbol_de(:,:,k),infor_symbol_de(:,:,k)] = MP_detection(K_sc_used,Rx_Sym_SC_Eq(:, :, k),2^Bit_Loading_SC(1),N_Symbol,abs(mean(conj(CRdata)./abs(CRdata).^2))^2*N_sigma,bit);
    end
 end

 
 %% 符号转二进制比特序列
Rx_bin_index = zeros(N_Symbol, N_Bit_Symbol_index);
Rx_bin_index = de2bi(mod(index_symbol_de,2^ N_Bit_Symbol_index), N_Bit_Symbol_index,'left-msb');%从左开始由二进制比特转化为符号
for i = 1:K_sc_used
    Rx_bin_info(:,Mod_Bit*i-Mod_Bit+1:Mod_Bit*i) = de2bi(infor_symbol_de(:,i), Mod_Bit,'left-msb'); 
end
Rx_bin = [Rx_bin_index Rx_bin_info];

output = Rx_bin;
 
end


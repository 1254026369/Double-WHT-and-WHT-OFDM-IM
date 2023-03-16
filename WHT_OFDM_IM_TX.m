function [output,output2,index_bit,info_bit] = WHT_OFDM_IM_TX(OFDM_IM_Para,bit)
%% OFDM_TX  OFDM_IM传输机

%% ---------- OFDM Parameters
N_Symbol = OFDM_IM_Para.N_Symbol;
N_frm = OFDM_IM_Para.N_frm;
Nb_fft = OFDM_IM_Para.Nb_fft;
CP = OFDM_IM_Para.CP;
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


%% 基带数据数据产生
P_data = bits_test(N_Bit_Symbol,N_Symbol,N_frm);
Date_bits = reshape(P_data,N_frm,[]);
save('Date_bits','Date_bits');

%% 组合法产生索引映射
index_all = Combin_Md(N_sc_used,K_sc_used);     % Combin_Md函数索引映射方式的实现    

%% OFDM IM调制(N_Symbol个信号并行同时调制)

Tx_Sym_FFT = zeros(N_Symbol, N_fft, N_frm);
for k = 1: N_frm
    %信息比特QAM调制映射
    Data_Bin = Date_bits(k, :);
    ss_mat = vec2mat(Data_Bin, N_Bit_Symbol);%将向量拆开转化为矩阵，N_Bit_Symbol列,最后补零
    ss_Mapping = zeros(N_Symbol, N_fft);   %N_fft列
    info_bit = ss_mat(:,N_Bit_Symbol_index+1:end);    %提取信息比特
    k_start = 1;
    
    %测试变量
    xxx = 1; 
    
    for k_SC = 1: N_fft
        Mod_Bit = Bit_Loading_SC_active(k_SC);
        if Mod_Bit == 0    % dropped subcarrier
            ss_Mapping(:, k_SC) = zeros(N_Symbol, 1);
        else               % adaptive modulation
            k_end = k_start + Mod_Bit - 1;
            bin_sc = info_bit(:, k_start: k_end);   
            dec_sc = bi2de(bin_sc, 'left-msb');%从左开始由二进制比特转化为符号
            ss_Mapping(:, k_SC) = qammod(dec_sc, 2^Mod_Bit, 'gray');%格雷编码
            k_start = k_end + 1;
            %测试语句
            tx_symbol(:,xxx) = dec_sc;
            xxx = xxx + 1;
            
        end  
    end
    index_bit = ss_mat(:,1:N_Bit_Symbol_index);    %提取索引比特
    index_sym = bi2de(index_bit, 'left-msb');%从左开始由二进制比特转化为符号
    tx_sym = zeros(N_Symbol,N_fft);     %N_sc_used载波，N_fft组的比特块
    for kk = 1: N_Symbol
            kk_index = index_sym(kk)+1;     %索引符号加一（matlab数组的索引为正整数）
            indices = index_all(kk_index,:)+1;      %索引符号映射到第indices 个子载波
            in_index = Find_index(indices,Bit_index);
            tx_sym(kk,in_index) = ss_Mapping(kk, find(Bit_Loading_SC_active));      %第kk组比特块选择第indices子载波进行传送该符号（比特）
    end
    Tx_Sym_FFT(:, :, k) = tx_sym; 
end

% 测试保存数据
save('index_bit','index_bit');
save('index_sym','index_sym');
save('tx_symbol','tx_symbol');
save('info_bit','info_bit');
save('Tx_Sym_FFT', 'Tx_Sym_FFT'); % save transmitted complex matrix

%% 功率归一化
% --------------Power-loading/Tx_Sym_FFT-----------------%
for k = 1: N_frm
    for k_SC = 1: N_fft
        Mod_Bit = Bit_Loading_SC(k_SC);
        switch Mod_Bit
            case 0
                P_mean = 1;
            case 2
                P_mean = 2;
            case 4
                P_mean = 10;     % normalized powers for different QAM modulation formats归一化功率
            case 5
                P_mean = 20;
            case 6
                P_mean = 42;
            case 7
                P_mean = 82;
            case 8
                P_mean = 170;
            otherwise
                error('Incorrect Subcarrier Modulation Format!')
        end 
        Tx_Sym_FFT(:, k_SC, k) = Tx_Sym_FFT(:, k_SC, k)/sqrt(P_mean); 
    end
end

%% WHT
Hada = hadamard(N_sc_used);
Tx_Sym_FFT_WHT =  Tx_Sym_FFT * Hada; 


%% IFFT + 循环前缀

if mod(N_fft * CP,N_fft) == 0
    OFDM_Sample = zeros(N_frm, N_Symbol * N_fft * (1 + CP));
else
    OFDM_Sample = zeros(N_frm, N_Symbol * N_fft + ceil(CP * N_fft) * N_Symbol);
end


for k = 1: N_frm
    FFT_Window = Tx_Sym_FFT(:, :, k);   
    OFDM_Window = ifft(FFT_Window, N_fft, 2);   % row wise FFT 每一行的FFT
    N_CP = ceil(CP * N_fft);                        % sample length of cyclic prefix 循环前缀采样长度
    CP_Samples = OFDM_Window(:, N_fft - N_CP + 1: N_fft);  % cyclic prefix samples 取出循环前缀
    OFDM_Window_CP = [CP_Samples, OFDM_Window];    % add CP
    OFDM_Signal = reshape(OFDM_Window_CP.', 1, numel(OFDM_Window_CP));     % change to a row vector 
    % numel()数组元素的数目
    % 10帧OFDM帧（6个）符号
    OFDM_Sample(k, :) = OFDM_Signal;    
end

for k = 1: N_frm
    FFT_Window_WHT = Tx_Sym_FFT_WHT(:, :, k);   
    OFDM_Window_WHT = ifft(FFT_Window_WHT, N_fft, 2);   % row wise FFT 每一行的FFT
    N_CP = ceil(CP * N_fft);                        % sample length of cyclic prefix 循环前缀采样长度
    CP_Samples_WHT = OFDM_Window_WHT(:, N_fft - N_CP + 1: N_fft);  % cyclic prefix samples 取出循环前缀
    OFDM_Window_CP_WHT = [CP_Samples_WHT, OFDM_Window_WHT];    % add CP
    OFDM_Signal_WHT = reshape(OFDM_Window_CP_WHT.', 1, numel(OFDM_Window_CP_WHT));     % change to a row vector 
    % numel()数组元素的数目
    % 10帧OFDM帧（6个）符号
    OFDM_Sample_WHT(k, :) = OFDM_Signal_WHT;    
end


% 限幅减小均峰比
% ------------OFDM Clipping/OFDM_LTE_Sample---------------%

if ClippingRatio ~= 0 % w/ clipping    
    for k = 1: N_frm               
        OFDM_Signal = OFDM_Sample_WHT(k, :);            
        OFDM_I_Signal = real(OFDM_Signal);  
        OFDM_I_Signal = OFDM_I_Signal - mean(OFDM_I_Signal); 
        
        %w/ clipping
        I_P_mean = mean(abs(OFDM_I_Signal).^2);
        I_V_max = sqrt(10^(ClippingRatio/10)*I_P_mean);    % Clipping stage, this defines the maximum amplitude values
        
%         I_V_max = ClippingRatio * max(abs(OFDM_I_Signal));
    
        I_Comp = (OFDM_I_Signal > I_V_max)*I_V_max + (OFDM_I_Signal < (-I_V_max))*(-I_V_max) + ...
             ((OFDM_I_Signal <= I_V_max) & (OFDM_I_Signal >= (-I_V_max))).*OFDM_I_Signal; 
                 

        OFDM_Q_Signal = imag(OFDM_Signal);  
        OFDM_Q_Signal = OFDM_Q_Signal - mean(OFDM_Q_Signal); %remove the DC component 除去直流分量
        
       % w/ clipping
        Q_P_mean = mean(abs(OFDM_Q_Signal).^2);
        Q_V_max = sqrt(10^(ClippingRatio/10)*Q_P_mean);    % Clipping stage, this defines the maximum amplitude values
        
%         Q_V_max = ClippingRatio * max(abs(OFDM_Q_Signal));
        
        Q_Comp = (OFDM_Q_Signal > Q_V_max)*Q_V_max + (OFDM_Q_Signal < (-Q_V_max))*(-Q_V_max) + ...
             ((OFDM_Q_Signal <= Q_V_max) & (OFDM_Q_Signal >= (-Q_V_max))).*OFDM_Q_Signal; 
      
        OFDM_Sample_WHT(k, :) = I_Comp + 1i*Q_Comp;
               
    end   
end  

%% 限幅减小均峰比
% if Clipping_Agg ~= 0
%     
%     Sig_Temp = OFDM_Sample;
%     
%      P_Agg_mean = mean(abs(Sig_Temp).^2,2);
%      V_max = sqrt(10^(Clipping_Agg/10)*P_Agg_mean);    % Clipping stage, which defines the maximum amplitude values
% 
%     Sig_Upsampled = (Sig_Temp > V_max).*V_max + (Sig_Temp < (-V_max)).*(-V_max) + ...
%                         ((Sig_Temp <= V_max) & (Sig_Temp >= (-V_max))).*Sig_Temp; 
% end
%% ------------DAC---------------%
N_bits = Quantizationbit; % 8-bits DAC/ADC
OFDM_I_Signal = real(OFDM_Signal);  
OFDM_I_Signal = OFDM_I_Signal - mean(OFDM_I_Signal);         
I_V_max = max(abs(OFDM_I_Signal));
OFDM_Q_Signal = imag(OFDM_Signal);  
OFDM_Q_Signal = OFDM_Q_Signal - mean(OFDM_Q_Signal); %remove the DC component 除去直流分量
Q_V_max = max(abs(OFDM_Q_Signal));
% Sig_Upsampled_r = real(Sig_Upsampled);
% Sig_Upsampled_i = imag(Sig_Upsampled);
dv_r = I_V_max*2/(2^N_bits);
dv_i = Q_V_max*2/(2^N_bits);
SigOut_DAC_r = (floor(OFDM_I_Signal./ dv_r) + 1/2) .* dv_r;
SigOut_DAC_i = (floor(OFDM_Q_Signal ./ dv_i ) + 1/2) .* dv_i ;
SigOut_DAC = SigOut_DAC_r+SigOut_DAC_i*j;
% d_V = 2 * V_max / (2^N_bits);
% SigOut_DAC = (floor(Sig_Upsampled ./ d_V) + 1/2) .* d_V;
output = SigOut_DAC;


OFDM_I_Signal_WHT = real(OFDM_Signal_WHT);  
OFDM_I_Signal_WHT = OFDM_I_Signal_WHT - mean(OFDM_I_Signal_WHT);         
I_V_max_WHT = max(abs(OFDM_I_Signal_WHT));
OFDM_Q_Signal_WHT = imag(OFDM_Signal_WHT);  
OFDM_Q_Signal_WHT = OFDM_Q_Signal_WHT - mean(OFDM_Q_Signal_WHT); %remove the DC component 除去直流分量
Q_V_max_WHT = max(abs(OFDM_Q_Signal_WHT));
% Sig_Upsampled_r = real(Sig_Upsampled);
% Sig_Upsampled_i = imag(Sig_Upsampled);
dv_r_WHT = I_V_max_WHT*2/(2^N_bits);
dv_i_WHT = Q_V_max_WHT*2/(2^N_bits);
SigOut_DAC_r_WHT = (floor(OFDM_I_Signal_WHT./ dv_r_WHT) + 1/2) .* dv_r_WHT;
SigOut_DAC_i_WHT = (floor(OFDM_Q_Signal_WHT ./ dv_i_WHT ) + 1/2) .* dv_i_WHT ;
SigOut_DAC_WHT = SigOut_DAC_r_WHT+SigOut_DAC_i_WHT*i;
% d_V = 2 * V_max / (2^N_bits);
% SigOut_DAC = (floor(Sig_Upsampled ./ d_V) + 1/2) .* d_V;
output2 = SigOut_DAC_WHT;

end

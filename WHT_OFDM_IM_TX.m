function [output,output2,index_bit,info_bit] = WHT_OFDM_IM_TX(OFDM_IM_Para,bit)
%% OFDM_TX  OFDM_IM�����

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


%% �����������ݲ���
P_data = bits_test(N_Bit_Symbol,N_Symbol,N_frm);
Date_bits = reshape(P_data,N_frm,[]);
save('Date_bits','Date_bits');

%% ��Ϸ���������ӳ��
index_all = Combin_Md(N_sc_used,K_sc_used);     % Combin_Md��������ӳ�䷽ʽ��ʵ��    

%% OFDM IM����(N_Symbol���źŲ���ͬʱ����)

Tx_Sym_FFT = zeros(N_Symbol, N_fft, N_frm);
for k = 1: N_frm
    %��Ϣ����QAM����ӳ��
    Data_Bin = Date_bits(k, :);
    ss_mat = vec2mat(Data_Bin, N_Bit_Symbol);%��������ת��Ϊ����N_Bit_Symbol��,�����
    ss_Mapping = zeros(N_Symbol, N_fft);   %N_fft��
    info_bit = ss_mat(:,N_Bit_Symbol_index+1:end);    %��ȡ��Ϣ����
    k_start = 1;
    
    %���Ա���
    xxx = 1; 
    
    for k_SC = 1: N_fft
        Mod_Bit = Bit_Loading_SC_active(k_SC);
        if Mod_Bit == 0    % dropped subcarrier
            ss_Mapping(:, k_SC) = zeros(N_Symbol, 1);
        else               % adaptive modulation
            k_end = k_start + Mod_Bit - 1;
            bin_sc = info_bit(:, k_start: k_end);   
            dec_sc = bi2de(bin_sc, 'left-msb');%����ʼ�ɶ����Ʊ���ת��Ϊ����
            ss_Mapping(:, k_SC) = qammod(dec_sc, 2^Mod_Bit, 'gray');%���ױ���
            k_start = k_end + 1;
            %�������
            tx_symbol(:,xxx) = dec_sc;
            xxx = xxx + 1;
            
        end  
    end
    index_bit = ss_mat(:,1:N_Bit_Symbol_index);    %��ȡ��������
    index_sym = bi2de(index_bit, 'left-msb');%����ʼ�ɶ����Ʊ���ת��Ϊ����
    tx_sym = zeros(N_Symbol,N_fft);     %N_sc_used�ز���N_fft��ı��ؿ�
    for kk = 1: N_Symbol
            kk_index = index_sym(kk)+1;     %�������ż�һ��matlab���������Ϊ��������
            indices = index_all(kk_index,:)+1;      %��������ӳ�䵽��indices �����ز�
            in_index = Find_index(indices,Bit_index);
            tx_sym(kk,in_index) = ss_Mapping(kk, find(Bit_Loading_SC_active));      %��kk����ؿ�ѡ���indices���ز����д��͸÷��ţ����أ�
    end
    Tx_Sym_FFT(:, :, k) = tx_sym; 
end

% ���Ա�������
save('index_bit','index_bit');
save('index_sym','index_sym');
save('tx_symbol','tx_symbol');
save('info_bit','info_bit');
save('Tx_Sym_FFT', 'Tx_Sym_FFT'); % save transmitted complex matrix

%% ���ʹ�һ��
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
                P_mean = 10;     % normalized powers for different QAM modulation formats��һ������
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


%% IFFT + ѭ��ǰ׺

if mod(N_fft * CP,N_fft) == 0
    OFDM_Sample = zeros(N_frm, N_Symbol * N_fft * (1 + CP));
else
    OFDM_Sample = zeros(N_frm, N_Symbol * N_fft + ceil(CP * N_fft) * N_Symbol);
end


for k = 1: N_frm
    FFT_Window = Tx_Sym_FFT(:, :, k);   
    OFDM_Window = ifft(FFT_Window, N_fft, 2);   % row wise FFT ÿһ�е�FFT
    N_CP = ceil(CP * N_fft);                        % sample length of cyclic prefix ѭ��ǰ׺��������
    CP_Samples = OFDM_Window(:, N_fft - N_CP + 1: N_fft);  % cyclic prefix samples ȡ��ѭ��ǰ׺
    OFDM_Window_CP = [CP_Samples, OFDM_Window];    % add CP
    OFDM_Signal = reshape(OFDM_Window_CP.', 1, numel(OFDM_Window_CP));     % change to a row vector 
    % numel()����Ԫ�ص���Ŀ
    % 10֡OFDM֡��6��������
    OFDM_Sample(k, :) = OFDM_Signal;    
end

for k = 1: N_frm
    FFT_Window_WHT = Tx_Sym_FFT_WHT(:, :, k);   
    OFDM_Window_WHT = ifft(FFT_Window_WHT, N_fft, 2);   % row wise FFT ÿһ�е�FFT
    N_CP = ceil(CP * N_fft);                        % sample length of cyclic prefix ѭ��ǰ׺��������
    CP_Samples_WHT = OFDM_Window_WHT(:, N_fft - N_CP + 1: N_fft);  % cyclic prefix samples ȡ��ѭ��ǰ׺
    OFDM_Window_CP_WHT = [CP_Samples_WHT, OFDM_Window_WHT];    % add CP
    OFDM_Signal_WHT = reshape(OFDM_Window_CP_WHT.', 1, numel(OFDM_Window_CP_WHT));     % change to a row vector 
    % numel()����Ԫ�ص���Ŀ
    % 10֡OFDM֡��6��������
    OFDM_Sample_WHT(k, :) = OFDM_Signal_WHT;    
end


% �޷���С�����
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
        OFDM_Q_Signal = OFDM_Q_Signal - mean(OFDM_Q_Signal); %remove the DC component ��ȥֱ������
        
       % w/ clipping
        Q_P_mean = mean(abs(OFDM_Q_Signal).^2);
        Q_V_max = sqrt(10^(ClippingRatio/10)*Q_P_mean);    % Clipping stage, this defines the maximum amplitude values
        
%         Q_V_max = ClippingRatio * max(abs(OFDM_Q_Signal));
        
        Q_Comp = (OFDM_Q_Signal > Q_V_max)*Q_V_max + (OFDM_Q_Signal < (-Q_V_max))*(-Q_V_max) + ...
             ((OFDM_Q_Signal <= Q_V_max) & (OFDM_Q_Signal >= (-Q_V_max))).*OFDM_Q_Signal; 
      
        OFDM_Sample_WHT(k, :) = I_Comp + 1i*Q_Comp;
               
    end   
end  

%% �޷���С�����
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
OFDM_Q_Signal = OFDM_Q_Signal - mean(OFDM_Q_Signal); %remove the DC component ��ȥֱ������
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
OFDM_Q_Signal_WHT = OFDM_Q_Signal_WHT - mean(OFDM_Q_Signal_WHT); %remove the DC component ��ȥֱ������
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

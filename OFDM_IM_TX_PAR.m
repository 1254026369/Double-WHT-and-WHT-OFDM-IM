clc;
clear;

%% ---------------- Control Parameters ---------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_Symbol = 10000;                         % Number of OFDM symbols to be transmitted   (OFDM_Para)需要传递的OFDM的符号个数
CPshift = 0;                            % CP offset for OFDM demodulation  (OFDM_Para)
N_frm = 1;          %帧数、frame
detection = 3;      %1:ML判决     2:LLR判决    3：MP

%% ---------------- OFDM IM Parameters ------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_QAM = 6;                          % 一个QAM符号所携带的比特个数
M = 2^M_QAM ; % bits per M-ary symbol M进制符号的比特数
BW_OFDM = 30.72e6;                   % LTE sampling rate LTE采样率
Nb_fft = 1;                          % oversampling ratio when generating the OFDM signal OFDM信号的过采样率
CP = 1/8;                            % length ratio of the cyclic prefix 循环前缀的长度比
ClippingRatio = 0;                   % PAPR coefficient: '0'- w/o clipping  峰均比系数：0代表不限幅
Clipping_Agg = 22; % Clipping_Agg is 14dB.限幅14dB

N_sc_used = 22;                      % data-carrying subcarrier number   数据携带载波个数
K_sc_used = 14;                  % data-carrying subcarrier number   数据携带激活子载波个数
N_fft = 2^nextpow2(N_sc_used);       % FFT/IFFT size  FFT/IFFT的大小
Quantizationbit = 11;

%按照惯例，nextpow2(0) 返回零。
%您可以使用 nextpow2 填充传递到 fft 的信号。
%当信号长度并非 2 次幂时，这样做可以加快 FFT 的运算速度。

%---- Bitloading
Bit_Loading_SC = zeros(1, N_fft);
Bit_Loading_SC_active = Bit_Loading_SC;
Bit_Loading_SC(1: N_sc_used) = M_QAM * ones(1, N_sc_used);
Bit_Loading_SC = circshift(Bit_Loading_SC, [0, floor(-N_sc_used/2)]);
Bit_Loading_SC_active(1: K_sc_used) = M_QAM * ones(1, K_sc_used);
Bit_Loading_SC_active= circshift(Bit_Loading_SC_active, [0, floor(-K_sc_used/2)]);
Bit_index = sort((1:N_fft) - (N_fft - N_sc_used),"descend");
Bit_index = circshift(Bit_index, [0, floor(-N_sc_used/2)]);

N_Bit_Symbol_info = sum(Bit_Loading_SC_active);         %number of transmitted information bits in each OFDM IM symbol
N_Bit_Symbol_index = floor(log2(nchoosek(N_sc_used,K_sc_used)));        %number of transmitted index bits in each OFDM IM symbol
N_Bit_Symbol = N_Bit_Symbol_info + N_Bit_Symbol_index;                % number of transmitted bits in each OFDM IM symbol每个OFDM符号传输的比特数量

OFDM_IM_Para = struct('N_Symbol', N_Symbol, 'CPshift', CPshift,'N_frm',N_frm, 'BW_OFDM', BW_OFDM,'detection',detection,...
                   'Nb_fft', Nb_fft, 'CP', CP, 'ClippingRatio', ClippingRatio, 'N_sc_used', N_sc_used,'K_sc_used',K_sc_used,...
                   'N_fft', N_fft, 'Bit_Loading_SC', Bit_Loading_SC,'Bit_Loading_SC_active', Bit_Loading_SC_active,'N_Bit_Symbol_info',N_Bit_Symbol_info,...
                   'N_Bit_Symbol_index',N_Bit_Symbol_index,'N_Bit_Symbol', N_Bit_Symbol,'Clipping_Agg',Clipping_Agg,'Bit_index',Bit_index,'Quantizationbit',Quantizationbit);
              
save('OFDM_IM_Para', 'OFDM_IM_Para'); 


%% ------------------ Tx Signal Generation -----------------%
[index_bit,index_sym,tx_symbol,info_bit,Tx_Signal_OS] = OFDM_IM_TX(OFDM_IM_Para);
save('Tx_Signal_OS','Tx_Signal_OS'); 
FFT_Tx_Signal_OS = fft(Tx_Signal_OS ,65536);
FFT_Tx_Signal_OS =fftshift(FFT_Tx_Signal_OS);
figure;
plot(-65535/2:65535/2,abs(FFT_Tx_Signal_OS));
DB=10*log10((abs(FFT_Tx_Signal_OS).^2)/0.001);
figure;
plot(-65535/2:65535/2,DB);

clc;
clear;

%% 载入数据
Rx_Sample = load('Tx_Signal_OS.mat');
load('OFDM_IM_Para.mat');
load('Date_bits.mat');

%% 接收机输出
[Rx_bin_index ,Rx_bin_info,index_symbol_de,infor_symbol_de] = OFDM_IM_RX(OFDM_IM_Para, Rx_Sample.Tx_Signal_OS,1);


function P_data = bits_test(N_Bit_Symbol,N_Symbol,Tx_SB_Index)
%% 基带数据数据产生
P_data=randi([0 1],1,N_Bit_Symbol*N_Symbol*Tx_SB_Index);
end


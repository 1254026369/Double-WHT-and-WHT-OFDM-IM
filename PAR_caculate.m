function PARdB = PAR_caculate(Tx_Signal_OS,N_fft)

Tx_Signal_OS = Tx_Signal_OS.*sqrt(N_fft);
PARdB = 10.*log10(max(abs(Tx_Signal_OS.^2),[],2)./(mean(abs(Tx_Signal_OS.^2),2)));

end


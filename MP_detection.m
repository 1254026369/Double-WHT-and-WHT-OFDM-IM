function [index_sym,info_sym_de] = N_MP_detection(K,y,M,N_Symbol,N_sigma,a)
%ML_DETECTION 最大功率估计
N = size(y,2);
index = zeros(N_Symbol,N);
index_sym = zeros(N_Symbol,1);
info_sym = zeros(N_Symbol,K);
info_sym_de = zeros(N_Symbol,K);

P_de = abs(y).^2;
[after_sort_P, index_P] = sort(P_de.','descend');
% for ii = 1:N_Symbol
%     FLAG = sum(after_sort_P(1:K,ii)>N_sigma*a);
%     if  FLAG~=K
%         for jj = 1:FLAG
%             index_P(jj + K - FLAG,ii) = index_P(jj + K,ii);
%         end
%     end        
% end
index = (N - index_P).';
index_P = sort(index_P(1:K,:)).';
index = sort(index(:,1:K).' ,'descend').';

    
for ii =1:N_Symbol
    info_sym(ii,:) = y(ii,index_P(ii,:));
end

info_sym_de = qamdemod(info_sym, M, 'gray');   
    
for ii = 1:N_Symbol
    KK = K;
    for jj =1:K
    index_sym(ii,1) = index_sym(ii,1) + binomial(index(ii,jj),KK);
    KK = KK - 1;
    end
end
end
function [index_sym,info_sym_de] = LLR_detection(K,y,h,M,N_Symbol,sym_test,N_sigma)
%LLR_DETECTION LLRÈíÅÐ¾ö
N = size(y,2);
lamda = zeros(N_Symbol,N);
index = zeros(N_Symbol,N);
index_sym = zeros(N_Symbol,1);
info_sym = zeros(N_Symbol,K);
info_sym_de = zeros(N_Symbol,K);

for ii = 1:N_Symbol
    re_LLR = repmat(y(ii,:),M,1);
    sym_de = repmat(sym_test.',1,N);
    exponential = exp((-abs((re_LLR-h.*sym_de)).^2)./N_sigma);
    lamda(ii,:) = log(K/(N-K)) + abs(y(ii,:)).^2/N_sigma + log(sum(exponential));
end

[after_sort_lamda, index_lamda] = sort(lamda.','descend');
index = (N - index_lamda).';
index_lamda = sort(index_lamda(1:K,:)).';
index = sort(index(:,1:K).' ,'descend').';

for ii =1:N_Symbol
    info_sym(ii,:) = y(ii,index_lamda(ii,:));
end

info_sym_de = qamdemod(info_sym, M, 'gray');   

for ii = 1:N_Symbol
    KK = K;
    for jj =1:K
    index_sym(ii,1) = index_sym(ii,1) + binomial(index(ii,jj),KK);
    KK = KK - 1;
    end
end


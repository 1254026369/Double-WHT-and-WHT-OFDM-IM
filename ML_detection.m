function [BB,MM] = ML_detection(K,N_Bit_Symbol_index,index_all,y,h,M,N_Symbol,sym_test)
%ML_DETECTION 最大似然估计
N = size(y,2);
for ii = 1:N_Symbol
    dis_id = zeros(1,2^N_Bit_Symbol_index);
    dis_sym = zeros(K,2^N_Bit_Symbol_index);    %K行2^p1列的矩阵（所有索引符号可能性的数列）
    
    for bb=1:2^N_Bit_Symbol_index
        id = index_all(bb,:)+1;   %所有索引符号中第bb个符号+1（目的是matlab的索引从1开始）
        id = N - id + 1;
        sym_k = zeros(K,1);   %K行1列
        sum = 0;
        for k1=1:K
            n1 = id(k1);    %第k1组符号映射
            disA = zeros(1,M);      %生成1行M列矩阵
            for k2=1:M
                % Power reallocation 功率再分配（对发送机输出功率）
                sym_m = sym_test(k2);
                disA(k2)=norm(y(ii,n1)-h*sym_m).^2;
            end
            [minA,iA] = min(disA);
            sum = sum + minA;   
            sym_k(k1) = sym_test(iA);     %ref_sym未将能量归一的
        end
        % dis_id(bb)=sum;
        dis_sym(:,bb)=sym_k; %%
        %%  New Low-Com ML detector
        sym_b = zeros(N,1);
        sym_b(id) = sym_k; %
        tmp = norm(y(ii,:).'-h*sym_b).^2;
        dis_id(bb)=tmp;
    end
    [~, I] = min(dis_id);
    BB(ii,:) = I - 1; 
    MM(ii,:) = dis_sym(:,I).';
end
end


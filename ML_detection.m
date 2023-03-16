function [BB,MM] = ML_detection(K,N_Bit_Symbol_index,index_all,y,h,M,N_Symbol,sym_test)
%ML_DETECTION �����Ȼ����
N = size(y,2);
for ii = 1:N_Symbol
    dis_id = zeros(1,2^N_Bit_Symbol_index);
    dis_sym = zeros(K,2^N_Bit_Symbol_index);    %K��2^p1�еľ��������������ſ����Ե����У�
    
    for bb=1:2^N_Bit_Symbol_index
        id = index_all(bb,:)+1;   %�������������е�bb������+1��Ŀ����matlab��������1��ʼ��
        id = N - id + 1;
        sym_k = zeros(K,1);   %K��1��
        sum = 0;
        for k1=1:K
            n1 = id(k1);    %��k1�����ӳ��
            disA = zeros(1,M);      %����1��M�о���
            for k2=1:M
                % Power reallocation �����ٷ��䣨�Է��ͻ�������ʣ�
                sym_m = sym_test(k2);
                disA(k2)=norm(y(ii,n1)-h*sym_m).^2;
            end
            [minA,iA] = min(disA);
            sum = sum + minA;   
            sym_k(k1) = sym_test(iA);     %ref_symδ��������һ��
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


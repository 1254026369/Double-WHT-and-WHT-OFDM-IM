function in_index = Find_index(indices,Bit_index)
%FIND_INDEX Ѱ�Ҷ�Ӧ�������ز���Ӧ����������
n = size(indices,2);
in_index = zeros(1,n);
for i = 1:n
    in_index(1,i) = find(Bit_index == indices(i));
end
end


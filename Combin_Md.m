
% Combinatorial Method 组合方式索引映射
% nn : number of subcarriers 子载波个数
% kk : number of active subcarriers 激活子载波个数
% output : all possible realization of active subcarrier indices 所有激活子载波的可能索引

function [output] = Combin_Md(nn,kk)

n = nn;
k = kk;

% all possible realization, binomial 所有可能实现的,二项式
N = binomial(n,k);

% For loop for all J sequence
JJ = zeros(N,k);
for ii = 1:N
    
    % Initialization 初始化
    n = nn; 
    k = kk;
    N = binomial(n,k);
    
    % Bit length of indices 索引的比特长度
    L = floor(log2(N));
    
    % J sequence J序列
    Z = binomial(n,k)-ii;
    tmp = zeros(1,k);
    
    % For loop to find combination sum of J(ii)
    for jj = 1:k
        
        ck = 0;
        C = 0;
        
        % While loop to find maximum of ck while C <=Z
        while C <= Z
            
            C = binomial(ck,k);
            
            ck = ck+1;
            
        end
        
        tmp(jj) = ck-2;
        % For next element
        Z = Z - binomial(tmp(jj),k);
        % k is second element of combination
        k = k - 1;
    end
    
    JJ(N-ii+1,:) = tmp;
    
end

output = JJ;
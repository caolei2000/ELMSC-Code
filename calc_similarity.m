function [sm] = calc_similarity(X1,X2,dim)
%CALC_SIMILARITY 计算两个视图数据间的相似度矩阵
%   此处显示详细说明

[~,n] = size(X1);
%% 使X1, X2维度保持一致
[~,score,~] = pca(X1');
X1 = score(:, 1:dim)';
[~,score,~] = pca(X2');
X2 = score(:, 1:dim)';
%% 使数据标准化(化为单位向量)
for i = 1:n
    X1(:,i) = X1(:,i) / norm(X1(:,i));
    X2(:,i) = X2(:,i) / norm(X2(:,i));
end
%% 计算余弦相似度
sm = zeros(n);
for i = 1:n
    for j = i:n
        sm(i,j) = dot(X1(:,i),X2(:,j));
    end    
end
sm = sm + sm' - diag(diag(sm));
sm = 1/2*sm + 0.5;  % MinMax归一化至[0,1]

end


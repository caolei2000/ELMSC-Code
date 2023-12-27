function [val] = obj_fun_val(E,nondiagJ,lambda,mu,Y1,Y2,Y3,X,P,L,Z,J,El,Ez)
%OBJ_FUN_VAL 计算BDLMSC的目标函数值
%   此处显示详细说明
E_l21 = mat_l21_norm(E);
J_l1 = sum(abs(nondiagJ),'all');
val = E_l21 + lambda*J_l1 + augmt_term(Y1,X-P*L-El,mu) + ...
    augmt_term(Y2,L-L*Z-Ez,mu) + augmt_term(Y3,J-Z,mu);

end

%% 辅助函数
function l21_norm = mat_l21_norm(A)
% 计算矩阵 A 的 L21 范数

% 计算每一列的 L2 范数
col_norms = sqrt(sum(A.^2, 1));
% 计算 L21 范数
l21_norm = sum(abs(col_norms));
end

function [val] = augmt_term(C,D,mu)
%AUGMT_TERM 计算增广拉格朗日函数的增广项和惩罚项
%   此处显示详细说明

val = mu/2*norm(D,"fro")^2 + trace(C'*D);

end


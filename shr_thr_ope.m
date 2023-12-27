function [y] = shr_thr_ope(x,alpha)
%SHR_THR_OPE 收缩阈值算子 shrinkage thresholding operator
%   此处显示详细说明

y = max(0,(abs(x)-alpha)) .* sign(x);

end


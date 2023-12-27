function P = update_P(Y1,mu,X,L,El)
% 闭式解更新P

Q = (1/mu*Y1+X-El)';
W = L*Q;
[U,~,V] = svd(W,'econ'); 
PT = U*V';
P = PT';

end

function [Z,Z_a,obj_val] = ELMSC(X,X_a,para)
% Enhanced Latent Multi-view Subspace Clustering, https://arxiv.org/abs/2312.14763v1
%Input:
%       X: multi-view data, cell, e.g., size(X{1})=d*n;
%       X_a: augmented data matrix;
%       para: algorithm parameters
%           K: the dimension of the latent representation;
%           lambda: balance parameter;
%Output:
%       Z: self-representation matrix, size=n*n;
%       Z_a: augmented self-representation matrix;
%       obj_val: objective function value
%
%   Written by Lei Cao 2023/12/27, leicao2000@gmail.com
%   
%% Initialization
V = size(X,2);  % Number of views
N = size(X{1},2);  % Number of samples
SD = 0;
SN = V*N;
for i = 1:V
    D(i) = size(X{i},1);  % Dimensions of each view
    SD = SD + D(i);
end
K = para.K;lambda = para.lambda;
rho = 1.2;mu = 1e-4;max_mu = 1e6;
P = zeros(SD,K);H_a = randn(K,SN);  
Z_a = zeros(SN,SN);J = zeros(SN,SN);
E1 = zeros(SD,SN);E2 = zeros(K,SN);
Y1 = zeros(SD,SN);Y2 = zeros(K,SN);Y3 = zeros(SN);
tol = 1e-2;
MaxIter = 10;
%% Start iterating with ADMM
for iter = 1:MaxIter
    % Update P
    P = update_P(Y1,mu,X_a,H_a,E1);
    % Update H_a
    A = mu*(P'*P);B = mu*(Z_a*Z_a'-Z_a-Z_a'+eye(SN));
    C = -(P'*Y1 + Y2*(Z_a'-eye(SN)) + mu*(P'*X_a+E2-P'*E1-E2*Z_a'));
    H_a = lyap(A,B,C);
    % Update Z_a
    Z_a = (H_a'*H_a+eye(SN)) \ ((J+H_a'*H_a-H_a'*E2)+(Y3+H_a'*Y2)/mu);
    Z_a = Z_a - diag(diag(Z_a));
    % Update E
    G = [X_a-P*H_a+Y1/mu;H_a-H_a*Z_a+Y2/mu];
    E = solve_l1l2(G,1/mu);
    E1 = E(1:SD,:);E2 = E(1+SD:SD+K,:);
    % Update J
    J = Z_a-Y3/mu;
    J1 = Z_a(1:N,1:N);
    for j = 2:V
        J1 = blkdiag(J1,J(j*N-N+1:j*N,j*N-N+1:j*N));
    end
    J2 = shr_thr_ope(J-J1,lambda/mu);
    J = J1+J2;
    % Update multipliers Y1,Y2,Y3
    Y1 = Y1 + mu*(X_a-P*H_a-E1);
    Y2 = Y2 + mu*(H_a-H_a*Z_a-E2);
    Y3 = Y3 + mu*(J-Z_a);
    mu = min(rho*mu,max_mu);
    % Record the objective function value
    obj_val(iter) = obj_fun_val(E,J2,lambda,mu,Y1,Y2,Y3,X_a,P,H_a,Z_a,J,E1,E2);
    % Check for convergence
    if (norm(X_a-P*H_a-E1,"inf")<tol && norm(H_a-H_a*Z_a-E2,"inf")<tol && norm(J-Z_a,"inf")<tol)
        break
    end
end
%% Calculate Z
Z = zeros(N);
for i = 1:V
    for j = 1:V
        Z = Z + Z_a(i*N-N+1:i*N,j*N-N+1:j*N);
    end
end
Z = Z - diag(diag(Z));
end


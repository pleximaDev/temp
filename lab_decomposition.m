clc
clear variables;

% Init
A = randn(4, 4);
b = randn(4, 1);

A = A * A'; % Making A positive-defined


%{
format long % short
format short
pi
c = pi - 5e-5
%}



%____________________________________________%
% Eigendecomposition of a matrix (spectral)
[V, D] = eig(A); %V - eigenvectors, D - eigenvalues
% [V2, D2, flag] = eigs(A) % differencies in positions of columns
Q = V;
L = D;
Q_inv = Q^(-1);

Ans = Q*Q_inv;
(Ans < 0.8);
Ans(Ans < 0.8) = 0;

disp(A)
B = Q * L * Q^(-1)
%____________________________________________%


%____________________________________________%
% Singular value decomposition
[U,S,V] = svd(A) 
% U - Unitary matrix with left-singular vectors;
% S - rectangular diagonal matrix with singular values;
% V - Unitary matrix with right-singular vectors in columns;
disp(A)
B = U * S * conj(V)'
%____________________________________________%


%____________________________________________%
% LU decomposition
[L,U] = lu(A)
disp(A)
B = L * U
%____________________________________________%


%____________________________________________%
% LUP decomposition
[L,U,P] = lu(A)
disp(A)
B = P' * L * U
%____________________________________________%


%____________________________________________%
% LL decomposition
L = chol(A,'lower')
U = chol(A, 'upper')
disp(A)
B = L * U
B = L * conj(L)'
B = conj(U)' * U
%____________________________________________%


%____________________________________________%
% LDL decomposition
[L,D] = ldl(A)
disp(A)
B = L * D * conj(L)'
%____________________________________________%


%____________________________________________%
% QR decomposition
[Q,R] = qr(A)
disp(A)
B = Q * R
%____________________________________________%



[L,U] = lu(A)
fprintf("_____________________")
[L,U] = my_lu(A, "Doolittle")
A
L*U

fprintf("_____________________\n")
[L,U] = my_lu(A, "Crout")
A
L*U
fprintf("_____________________\n")



[my, p] = my_chol(A, "Cholesky_Banachiewicz")
[standard] = chol(A, 'lower')
asymmentric_matrix = [1.4, 1, 2; 1, 0.9, 1; 1, 1, 1.4]
%%% Cholesky for asymmetric matrix %%%
[L_asymm, ok] = my_chol(asymmentric_matrix, "Cholesky_Banachiewicz")
B = L_asymm * conj(L_asymm)'
%%%               //               %%%
return

[Q_std, R_std] = qr(A)
[Q_classical, R] = my_qr(A, "Classical Gram_Schmidt")
[Q_modified, R] = my_qr(A, "Modified Gram_Schmidt")
[Q_Householder, R] = my_qr(A, "Householder")
[Q_Givens, R] = my_qr(A, "Givens")


function [L,U] = my_lu(x, method)
switch method
    case "Doolittle"
        [m, n] = size(x);
        L = zeros(m, n);
        U = zeros(m, n);
        sumL = 0;
        sumU = 0;
        for i = 1 : 1 : n
            for j = i : 1 : n
                for k = 1 : 1 : i-1
                    sumU = sumU + L(i, k) * U(k, j);
                    sumL = sumL + L(j, k) * U(k, i);
                end
                U(i,j) = x(i,j) - sumU;
                sumU = 0;
                L(i,i) = 1;
                if j ~= n
                    j = j + 1;
                end
                L(j,i) = (1/U(i,i)) * (x(j, i) - sumL);
                sumL = 0;
                j = j - 1;
            end    
        end
        
    case "Crout"
        [m, n] = size(x);
        L = zeros(m, n);
        U = zeros(m, n);
        sumL = 0;
        sumU = 0;
        for i = 1 : 1 : n
            for j = i : 1 : n
                for k = 1 : 1 : i-1
                    sumL = sumL + L(j, k) * U(k, i);
                    sumU = sumU + L(i, k) * U(k, j);
                end
                L(j,i) = x(j,i) - sumU;
                sumU = 0;
                U(i,i) = 1;
                if j ~= n
                    j = j + 1;
                end
                U(i,j) = (1/L(i,i)) * (x(i, j) - sumL);
                sumL = 0;
                j = j - 1;
            end    
        end
        
    otherwise 
        fprintf("Error occured while entering method's name.");
end
end

function [L, p] = my_chol(x, method)
switch method
    case "Cholesky_Banachiewicz"
        p = true;
        if all(eig(x) <= 1e-9)
            p = false;
            L = zeros(m, n);
            return
        end
        [m, n] = size(x);
        L = zeros(m, n);
        sumij = 0;
        sumii = 0;
        j = 1;
        k = 1;
        L(1, 1) = (x(1, 1))^(1/2);
        for i = 2 : 1 : m
            j = 1;
            while j <= i
               while k < j
                   sumij = sumij + L(i, k) * L(j, k);
                   k = k + 1;
               end
               k = 1;
               while k < j
                   sumii = sumii + L(i,k) * L(i,k);
                   k = k + 1;
               end
               k = 1;
               
               L(i,j) = (1/L(j,j))*(x(i,j) - sumij);
               L(i,i) = (x(i,i) - sumii)^(1/2);
               sumij = 0;
               sumii = 0;
               j = j + 1;
            end
        end
    otherwise
        fprintf("Error occured while entering method's name.")
end
end

function [Q,R] = my_qr(x, method)
switch method
    case "Classical Gram_Schmidt"
        [m,n] = size(x);
        Q = zeros(m,n);
        sum_proj = zeros(m,1);

        for k = 1 : 1 : m
            for j = 1 : 1 : k-1
                sum_proj = sum_proj + ((x(:,k)' * Q(:,j))/(Q(:,j)'*Q(:,j))).*Q(:,j);
            end
            Q(:,k) = x(:, k) - sum_proj;
            Q(:,k) = -Q(:,k)/norm(Q(:,k));
            sum_proj = 0;
        end
        R = Q' * x;
        %___________Unnecessary___________%
        zero_logic = triu(ones(m,n));
        logic = (zero_logic == 0);
        R(logic) = 0;
        %___________Unnecessary___________%
        
    case "Modified Gram_Schmidt"
        [m,n] = size(x);
        Q = zeros(m, n);
        b = zeros(m, n);

        b(:, 1) = x(:, 1)/norm(x(:, 1));
        for j = 2 : 1 : m
            b_prev = x(:, j);
            for i = 2 : 1 : j-1
                proj = ((b_prev' * b(:, i - 1))/(b(:, i - 1)'*b(:, i - 1))).*b(:, i - 1);
                b_curr = b_prev - proj;
                b_prev = b_curr;
            end
        % % % %     
        b(:, j) = b_prev - ((b_prev' * b(:, j - 1))/(b(:, j - 1)'*b(:, j - 1))).*b(:, j - 1);
        b(:, j) = b(:, j)/norm(b(:, j));
        end
        Q = -b;
        R = Q' * x;
        %___________Unnecessary___________%
        zero_logic = triu(ones(m,n));
        logic = (zero_logic == 0);
        R(logic) = 0;
        %___________Unnecessary___________%
    case "Householder"
        % Reflections
        A = x;
        [m, n] = size(A);
        Q = eye(m, n);
        e = zeros(m, 1);
        

        for k = 1 : 1 : n-1
            e(k, 1) = 1;
            x = A(:, k);
            for i = 1 : 1 : k-1
                x(i, 1) = 0;
            end
            u = x - norm(x) * e;
            P = eye(m, n) - 2 * (u * u')/(norm(u))^2;
            Q = Q * P;
            A = P * A;
            R = A;
            e(e>1e-12) = 0;
        end
        % <changing signs>
        R = -R;
        Q = -Q;
        % <\changing signs>
        %___________Unnecessary___________%
        zero_logic = triu(ones(m,n));
        logic = (zero_logic == 0);
        R(logic) = 0;
        %___________Unnecessary___________%
        
    case "Givens"
        % Rotations
        [m,n] = size(x);
        A = x;
        Q = eye(m, n);
        R = A;
        G = eye(n, n);
        for i = 1 : 1 : n - 1
            for j = i + 1 : 1 : n
                G = eye(n, n);
%                 fprintf("ij(%d, %d)\n", i, j);
                G(i, i) = A(i, i)/sqrt(A(i, i)^2 + A(j, i)^2);
                G(j, j) = G(i,i);
                G(i, j) = A(j, i)/sqrt(A(i, i)^2 + A(j, i)^2);
                G(j, i) = -A(j, i)/sqrt(A(i, i)^2 + A(j, i)^2);
                Q = Q * G';
                A = G * R;
                R = A;
            end
        end
        % <changing signs>
        R = -R;
        Q = -Q;
        % <\changing signs>
        %___________Unnecessary___________%
        zero_logic = triu(ones(m,n));
        logic = (zero_logic == 0);
        R(logic) = 0;
        %___________Unnecessary___________%
    otherwise 
        fprintf("Error occured while entering method's name.");
end
end







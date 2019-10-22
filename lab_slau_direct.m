clc
clear variables;
close all force;

%=========================================================================%
% #2

A = [1, -0.2589, -0.3093; -0.2589, 1, -0.2705; -0.3093, -0.2705, 1];
b = ones(3, 1);
x = [2.2873; 2.2162; 2.3068];


[x, ok] = lab_slau_gauss(A, b);

[x, ok] = lab_slau_gauss_jordan(A, b);

[x, ok] = lab_slau_minv(A, b);

[x, ok] = lab_slau_Cramer(A, b);

[x, ok] = lab_slau_chol(A, b);

%=========================================================================%
% #3

K = 4;
b = randn(K, 1);
I = eye(K);

%/\/\/\ A is symmetric positive definite diagonally dominant matrix  /\/\/\%
A0 = randn(K, K);
A0 = tril(A0);

A = A0 * A0' + 5 * K * I;
D{1} = {A, b, "Symmetric positive definite diagonally dominant matrix"};
%/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\%



%/\/\/\ A is symmetric negative definite diagonally dominant matrix  /\/\/\%
A0 = randn(K, K);
A0 = tril(A0);
A = A0 * A0' - 5 * K * I;
D{2} = {A, b, "Symmetric negative definite diagonally dominant matrix"};
%/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\%



%/\/\/\/\/\/\/\/\/\/\ A is non-symmetric random matrix /\/\/\/\/\/\/\/\/\/\%
A = randn(K, K);
D{3} = {A, b, "Non-symmetric random matrix"};
%/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\%



%/\/\/\/\ A is sparse negative definite diagonally dominant matrix /\/\/\/\%
%(singular, A^(-1) doesn't exist, detA = 0) n == m^2
A = gallery('neumann', K);
D{4} = {A, b, "Sparse negative definite diagonally dominant matrix"};
%/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\%




%===============================Saving data===============================%
save('lab_slau_data.mat', 'D', '-v7');
clear variables;
load('lab_slau_data.mat');
%=========================================================================%





%=========================================================================%
% #4

%################################################%
%@ X(4, 5, 4) - matrix of solutions             %@
%@ Where columns are solutions for each method  %@
%@ (:, :, i) are solutions for each matrix      %@
%@______________________________________________%@
%@                                              %@
%@ T(4, 5) - matrix of method's execution time  %@
%@ columns - methods                            %@
%@ row - matrices                               %@
%################################################%

for i = 1 : 1 : 4
    A = D{i}{1};
    b = D{i}{2};
    
    tic
    [x, ok] = lab_slau_gauss(A, b);
    t = toc
    T(i, 1) = t * ok;
    X(:, 1, i) = x;
    
    tic
    [x, ok] = lab_slau_gauss_jordan(A, b);
    t = toc
    T(i, 2) = t * ok;
    X(:, 2, i) = x;
    
    tic
    [x, ok] = lab_slau_minv(A, b);
    t = toc
    T(i, 3) = t * ok;
    X(:, 3, i) = x;
    
    tic
    [x, ok] = lab_slau_Cramer(A, b);
    t = toc
    T(i, 4) = t * ok;
    X(:, 4, i) = x;
    
    tic
    [x, ok] = lab_slau_chol(A, b);
    t = toc
    T(i, 5) = t * ok;
    X(:, 5, i) = x;
end

%=========================================================================%


N = 49
T = zeros(4, 5);

for i = 1 : 1 : 4
    A = D{i}{1};
    b = D{i}{2};
    
%     Vector_tmp = 0;
    for j = 1 : 1 : N
        tic
        [x, ok] = lab_slau_gauss(A, b);
        if ~ok
            break
        end
        Vector_tmp(j) = toc;
    end
    T(i, 1) = mean(Vector_tmp * ok);
    
    
%     Vector_tmp = 0;
    for j = 1 : 1 : N
        tic
        [x, ok] = lab_slau_gauss_jordan(A, b);
        if ~ok
            break
        end
        Vector_tmp(j) = toc;
    end
    T(i, 2) = mean(Vector_tmp * ok);
    
    
%     Vector_tmp = 0;
    for j = 1 : 1 : N
        tic
        [x, ok] = lab_slau_minv(A, b);
        if ~ok
            break
        end
        Vector_tmp(j) = toc;
    end
    T(i, 3) = mean(Vector_tmp * ok);
    
%     Vector_tmp = 0;
    for j = 1 : 1 : N
        tic
        [x, ok] = lab_slau_Cramer(A, b);
        if ~ok
            break
        end
        Vector_tmp(j) = toc;
    end
    T(i, 4) = mean(Vector_tmp * ok);
    
%     Vector_tmp = 0;
    for j = 1 : 1 : N
        tic
        [x, ok] = lab_slau_chol(A, b);
        if ~ok
            break
        end
        Vector_tmp(j) = toc;
    end
    T(i, 5) = mean(Vector_tmp * ok);
end



for i = 1 : 1 : 4
    Gauss = X(:, 1, i);
    Gauss_Jordan = X(:, 2, i);
    Cramer = X(:, 3, i);
    Invertible_matrix = X(:, 4, i);
    Cholesky = X(:, 5, i);
    fprintf("%s\n\n", D{i}{3});
    fprintf("A");
    disp(D{i}{1});
    fprintf("b");
    disp(D{i}{2});
    table(Gauss, Gauss_Jordan, Invertible_matrix, Cholesky);
end
% methods = ('Gauss','Gauss Jordan','Cramer','Invertible matrix','Cholesky');


clf;
subplot(1,2,1);
bar(T')
title('Direct methods')
ylabel('Time, ms');
xlabel('Methods');
ax = gca;
ax.XTickLabel={'Gauss','Gauss-Jordan',' Cramer',' Inverse ',' Cholesky'};
grid on
grid minor
legend({'A>0, Symmetric','A<0, Symmetric','A non-symmetric randn','A<0, Sparse'},'location','northeastoutside');


T(4, :) = []


subplot(1,2,2);
bar(T')
title('Direct methods')
ylabel('Time, ms');
xlabel('Methods');
ax = gca;
ax.XTickLabel={'Gauss','Gauss-Jordan',' Cramer',' Inverse ',' Cholesky'};
grid on
grid minor
legend({'A>0, Symmetric','A<0, Symmetric','A non-semmytric randn'},'location','northeastoutside');



%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
function [x, ok] = lab_slau_gauss(A, b)
ok = false;
[m, n] = size(A);
x = zeros(m, 1);

if n == m
    ok = true;
end

if ok 
    for i = 1 : 1 : n
        for k = i + 1 : 1 : m
            coefficient = (A(k, i)/A(i, i));
            A(k, :) = A(k, :) - A(i, :) * coefficient;
            b(k) = b(k) - b(i) * coefficient;
        end
    end
    for i = n : -1 : 1
        sum = 0;
        for k = i + 1 : 1 : n
            sum = sum + A(i, k) * x(k);
        end
        x(i) = (b(i) - sum)/A(i, i);
    end
end
end

function [x, ok] = lab_slau_gauss_jordan(A, b)
ok = false;
[m, n] = size(A);
x = zeros(m, 1);
C = horzcat(A, b);

if n == m
    ok = true;
end

if ok
    for i = 1 : 1 : n
        for k = i + 1 : 1 : n
            C(k, :) = C(k, :) - C(i, :) * (C(k, i)/C(i, i));
        end
    end
    for i = n : -1 : 1
        for k = i - 1 : -1 : 1
            C(k, :) = C(k, :) - C(i, :) * (C(k, i)/C(i, i));
        end
    end
    for t = 1 : 1 : m
        x(t) = C(t, m + 1)/C(t, t);
    end
    
end
end

function [x, ok] = lab_slau_Cramer(A, b)
[m,n] = size(A);
ok = false;
x = zeros(m, 1);
eps = 1e-16;

if m == n
    if det(A) < eps || det(A) > eps
        ok = true;
    end
end

if ok
    for i = 1 : 1 : n
        M = A;
        M(:, i) = b;
        x(i, 1) = (1/det(A)) * det(M);
    end
    
end
end

function [x, ok] = lab_slau_minv(A, b)
[m,n] = size(A);
ok = false;
x = zeros(m, 1);

if (m == n) && (-1e-16 > det(A) || 1e-16 < det(A))
    ok = true;
end

if ok
    for i = 1 : 1 : m
        for j = 1 : 1 : n
            M = A;
            M(i, :) = [];
            M(:, j) = [];
            C(i, j) = (-1)^(i + j) * det(M); 
        end
        T = (1 / det(A)) * C;
        x = T * b;
    end
end
end

function [x, ok] = lab_slau_chol(A, b)
ok = false;
[m, n] = size(A);
L = zeros(m, n);
x = zeros(m, 1);
sumij = 0;
sumii = 0;

if (all(eigs(A, m) > 1e-12))
    ok = true;
end

if ok 
    L(1, 1) = (A(1, 1))^(1/2);
    for i = 2 : 1 : m
        for j = 1 : 1 : i
           for k = 1 : 1 : j - 1
               sumij = sumij + L(i, k) * L(j, k);
           end

           L(i,j) = (1/L(j,j))*(A(i,j) - sumij);
           sumij = 0;
           sumii = 0;
        end
        for k = 1 : 1 : j - 1
               sumii = sumii + L(i,k) * L(i,k);
        end
        L(i,i) = (A(i,i) - sumii)^(1/2);
    end
    
    
    for i = 1 : 1 : n
        sum = 0;
        for k = 1 : 1 : i - 1
            sum = sum + L(i, k) * y(k);
        end
        y(i) = (b(i) - sum)/L(i, i);
    end
    for i = n : -1 : 1
        sum = 0;
        for k = i + 1 : 1 : n
            sum = sum + L(k, i) * x(k);
        end
        x(i) = (y(i) - sum)/L(i, i);
    end
    
%     y = L^(-1) * b;
%     x = (L^(-1))' * y;
end
end



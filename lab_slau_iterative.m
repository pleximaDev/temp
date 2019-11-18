clc
clear variables;
format short


%=========================================================================%
%=================================task #2=================================%
A = [1, -0.2589, -0.3093; -0.2589, 1, -0.2705; -0.3093, -0.2705, 1];
b = ones(3, 1);
x = [2.2873; 2.2162; 2.3068];
x0 = zeros(size(A, 1), 1);
eps = 1e-10;
Kmax = 1000;



[x, ok, k] = lab_slau_jacobi(A, b, x0, eps, Kmax);

[x, ok, k] = lab_slau_gauss_seidel(A, b, x0, eps, Kmax);

[x, ok, k, omega] = lab_slau_sor(A, b, x0, eps, Kmax);

[x, ok, k] = lab_slau_mcg(A, b, x0, eps, Kmax);

[x, ok, k] = lab_slau_mbcg(A, b, x0, eps, Kmax);

[x, ok, k] = lab_slau_smbcg(A, b, x0, eps, Kmax);
%=========================================================================%
%=========================================================================%



%=========================================================================%
%==============================obtaining data=============================%
K = 16;
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
x0 = zeros(K, 1);
D{5} = {x0, eps, Kmax, "Some required vars (x0, eps, Kmax)"};

%=========================================================================%
%=========================================================================%


%=========================================================================%
%===============================Saving data===============================%
save('lab_slau_data.mat', 'D', '-v7');
clear variables;
load('lab_slau_data.mat');
%=========================================================================%
%=========================================================================%




%=========================================================================%
%=============================time estimation=============================%

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

x0 = D{5}{1};
eps = D{5}{2};
Kmax = D{5}{3};
for i = 1 : 1 : 4
    A = D{i}{1};
    b = D{i}{2};
    
    tic
    [x, ok, k] = lab_slau_jacobi(A, b, x0, eps, Kmax);
    t = toc;
    T(i, 1) = t * ok;
    X(:, 1, i) = x;
    K(1, 1, i) = k;
    
    tic
    [x, ok, k] = lab_slau_gauss_seidel(A, b, x0, eps, Kmax);
    t = toc;
    T(i, 2) = t * ok;
    X(:, 2, i) = x;
    K(1, 2, i) = k;
   
    tic
    [x, ok, k] = lab_slau_sor(A, b, x0, eps, Kmax);
    t = toc;
    T(i, 3) = t * ok;
    X(:, 3, i) = x;
    K(1, 3, i) = k;
    
    tic
    [x, ok, k] = lab_slau_mcg(A, b, x0, eps, Kmax);
    t = toc;
    T(i, 4) = t * ok;
    X(:, 4, i) = x;
    K(1, 4, i) = k;
    
    tic
    [x, ok, k] = lab_slau_mbcg(A, b, x0, eps, Kmax);
    t = toc;
    T(i, 5) = t * ok;
    X(:, 5, i) = x;
    K(1, 5, i) = k;
    
    tic
    [x, ok, k] = lab_slau_smbcg(A, b, x0, eps, Kmax);
    t = toc;
    T(i, 6) = t * ok;
    X(:, 6, i) = x;
    K(1, 6, i) = k;
end



N = 1000;
Vector_tmp = 0;
for i = 1 : 1 : 4
    A = D{i}{1};
    b = D{i}{2};

    for j = 1 : 1 : N
        tic
        [x, ok, k] = lab_slau_jacobi(A, b, x0, eps, Kmax);
        if ~ok
            break
        end
        Vector_tmp(j) = toc;
    end
    T(i, 1) = mean(Vector_tmp * ok);
    

    for j = 1 : 1 : N
        tic
        [x, ok, k] = lab_slau_gauss_seidel(A, b, x0, eps, Kmax);
        if ~ok
            break
        end
        Vector_tmp(j) = toc;
    end
    T(i, 2) = mean(Vector_tmp * ok);
    

    for j = 1 : 1 : N
        tic
        [x, ok, k] = lab_slau_sor(A, b, x0, eps, Kmax);
        if ~ok
            break
        end
        Vector_tmp(j) = toc;
    end
    T(i, 3) = mean(Vector_tmp * ok);
    

    for j = 1 : 1 : N
        tic
        [x, ok, k] = lab_slau_mcg(A, b, x0, eps, Kmax);
        if ~ok
            break
        end
        Vector_tmp(j) = toc;
    end
    T(i, 4) = mean(Vector_tmp * ok);
    

    for j = 1 : 1 : N
        tic
        [x, ok, k] = lab_slau_mbcg(A, b, x0, eps, Kmax);
        if ~ok
            break
        end
        Vector_tmp(j) = toc;
    end
    T(i, 5) = mean(Vector_tmp * ok);
    
    for j = 1 : 1 : N
        tic
        [x, ok, k] = lab_slau_smbcg(A, b, x0, eps, Kmax);
        if ~ok
            break
        end
        Vector_tmp(j) = toc;
    end
    T(i, 6) = mean(Vector_tmp * ok);
end

%=========================================================================%
%=========================================================================%


for i = 1 : 1 : 4
    fprintf("=========================================================================\n");
    Jacobi = X(:, 1, i);
    Gauss_Seidel = X(:, 2, i);
    Successive_over_relaxation = X(:, 3, i);
    Conjugate_gradient = X(:, 4, i);
    Biconjugate_gradient = X(:, 5, i);
    Stabilized_biconjugate_gradient = X(:, 6, i);
    fprintf("%s\n\n", D{i}{3});
    fprintf("A");
    disp(D{i}{1});
    fprintf("b");
    disp(D{i}{2});
    table(Jacobi, Gauss_Seidel, Successive_over_relaxation, ...
    Conjugate_gradient, Biconjugate_gradient, Stabilized_biconjugate_gradient)
end
% Jacobi, Gauss_Seidel, Successive_over_relaxation, Conjugate_gradient, Biconjugate_gradient, Stabilized_biconjugate_gradient



figure('Name', 'Measured Data', 'NumberTitle', 'off');
clf;
subplot(2,1,1);
bar(T')
title('Iterative methods')
ylabel('Time, ms');
xlabel('Methods');
ax = gca;

ax.XTickLabel={'Jacobi', 'Gauss-Seidel', 'Successive over-relaxation', ...
'Conjugate gradient', 'Biconjugate gradient', 'Stabilized biconjugate gradient'};
grid on
grid minor
legend({'A>0, Symmetric','A<0, Symmetric','A non-symmetric randn','A<0, Sparse'},'location','northeastoutside');


%--------------------
%-----iterations-----
subplot(2, 1, 2);
bar(squeeze(K))
grid on
grid minor
ylabel('Iterations, k');
xlabel('Methods');
ax = gca;
ax.XTickLabel={'Jacobi', 'Gauss-Seidel', 'Successive over-relaxation', ...
'Conjugate gradient', 'Biconjugate gradient', 'Stabilized biconjugate gradient'};
legend({'A>0, Symmetric','A<0, Symmetric','A non-symmetric randn','A<0, Sparse'},'location','northeastoutside');
%--------------------
%--------------------






function [x, ok, k] = lab_slau_jacobi(A, b, x0, eps, Kmax)
if issparse(A)
    A = full(A);
end
[m, n] = size(A);
D = diag(diag(A));
L = tril(A, -1);
U = triu(A, 1);
k = 0;
x_k = x0;
x_kplus1 = b;

x = x0;
B = -D^(-1) * (L + U);
g = D^(-1) * b;

% if (eig(full(B)) < 1); ok = true; else; ok = false; end
ok = (all(eig(B) < 1));

if ok 
    while(k < Kmax) && (norm(x_kplus1 - x_k) > eps) %%%%%% FFFFFF
        x_k = x_kplus1;
        for i = 1 : 1 : m
            sum = 0;
            for j = 1 : 1 : n
                if j ~= i
                    sum = sum + A(i, j) * x_k(j);
                end
            end
            x_kplus1(i) = 1/A(i, i) * (b(i) - sum);
        end
        k = k + 1;
    end
end

x = x_k;

end

function [x, ok, k] = lab_slau_gauss_seidel(A, b, x0, eps, Kmax)
if issparse(A)
    A = full(A);
end
[m, n] = size(A);
D = diag(diag(A));
L = tril(A, -1);
U = triu(A, 1);
x = x0;
k = 0;

P = -(L + D)^(-1) * U;
ok = norm(P) < 1;

if ok
    for i = 1 : 1 : m
        for j = 1 : 1 : n
            if j ~= i
                c(i, j) = (-A(i, j)/A(i, i)) * (j ~= i) + 0 * (j == i);%%%FFFFF
            end
        end
    end
    
    while (k < Kmax) && (norm(A * x - b) > eps)
        for i = 1 : 1 : m
            d(i) = b(i)/A(i, i);
            sum1 = 0;
            for j = 1 : 1 : i - 1
                sum1 = sum1 + c(i, j) * x(j);
            end
            sum2 = 0;
            for j = i + 1 : 1 : n
                sum2 = sum2 + c(i, j) * x(j);
            end
            x(i) = sum1 + sum2 + d(i);
        end
        k = k + 1;
    end
end
end


function [x, ok, k, omega] = lab_slau_sor(A, b, x0, eps, Kmax)
% Successive over-relaxation
if issparse(A)
    A = full(A);
end
[m, n] = size(A);
D = diag(diag(A));
L = tril(A, -1);
U = triu(A, 1);
x = x0;
k = 0;
T = eye(m, n) - D^(-1) * A;
ro = max(abs(eig(T)));
omega = 1 - (ro/(1 + sqrt(1 - ro^2)))^2;

P = -(L + D)^(-1) * U;
ok = norm(P) < 1;

if ok
    while (k < Kmax) && (norm(A * x - b) > eps)
        for i = 1 : 1 : m
            sum1 = 0;
            for j = 1 : 1 : i - 1
                sum1 = sum1 + A(i, j) * x(j);
            end
            sum2 = 0;
            for j = i + 1 : 1 : n
                sum2 = sum2 + A(i, j) * x(j);
            end
            x(i) = (1 - omega) * x(i) + omega/A(i, i) * (b(i) - sum1 - sum2);
        end
        k = k + 1;
    end
end
end

function [x, ok, k] = lab_slau_mcg(A, b, x0, eps, Kmax)
% method of conjugate gradient
if issparse(A)
    A = full(A);
end
x = x0;
r = b - A * x;
z = r;
k = 0;

ok = isreal(A) * issymmetric(A) * (all(eig(A) > 0));

if ok
    while(k < Kmax) &&  (norm(A * x - b) > eps)
        r_prev = r;
        x_prev = x;
        z_prev = z;
        alpha = dot(r_prev, r_prev)/dot(A * z_prev, z_prev);
        x = x_prev + alpha * z_prev;
        r = r_prev - alpha * A * z_prev;
        betta = dot(r, r)/(dot(r_prev, r_prev));
        z = r + betta * z_prev;
        k = k + 1;
    end
end
end

function [x, ok, k] = lab_slau_mbcg(A, b, x0, eps, Kmax)
if issparse(A)
    A = full(A);
end
% method of biconjugate gradient
x = x0;
r = b - A * x0;
p = r;
z = r;
s = r;

k = 0;
ok = isreal(A) * (all(eig(A) > 0));

if ok
    while(k < Kmax) &&  (norm(A * x - b) > eps)
        p_prev = p;
        r_prev = r;
        x_prev = x;
        z_prev = z;
        s_prev = s;
        alpha = dot(p_prev, r_prev)/dot(s_prev, A * z_prev);
        x = x_prev + alpha * z_prev;
        r = r_prev - alpha * A * z_prev;
        p = p_prev - alpha * A' * s_prev;
        betta = dot(p, r)/dot(p_prev, r_prev);
        z = r + betta * z_prev;
        s = p + betta * s_prev;
        k = k + 1;
    end
end
end
    

function [x, ok, k] = lab_slau_smbcg(A, b, x0, eps, Kmax)
if issparse(A)
    A = full(A);
end
% Stabilized method of biconjugate gradient
my_dot_cmplx = @(q, d) conj(q') * d;
my_dot = @(q, d) q' * d;
k = 0;

x = x0;
r = b - A * x;
r_tilda = r;
ro = 1;
alpha = 1; 
omega = 1;
v = 0;
p = 0;

ok = isreal(A) * (all(eig(A) > 0));

if ok
    while(k < Kmax) &&  (norm(A * x - b) > eps)
        r_prev = r;
        ro_prev = ro;
        p_prev = p;
        alpha_prev = alpha;
        omega_prev = omega;
        v_prev = v;
        x_prev = x;
        
        ro = my_dot(r_tilda, r_prev);
        betta = (ro/ro_prev)*(alpha_prev/omega_prev);
        p = r_prev + betta * (p_prev - omega_prev * v_prev);
        v = A * p;
        alpha = ro/my_dot(r_tilda, v);
        s = r_prev - alpha * v;
        t = A * s;
        omega = my_dot_cmplx(t, s)/my_dot_cmplx(t, t);
        x = x_prev + alpha * p + omega * s;
        r = s - omega * t;
        k = k + 1;
    end
end
end


    
    
    
    

function [x, ok, k]=my_Jacobi(A, b, x0, kmax, epsilon)
% итерационный вид
% xk+1 = B * xk + g
% для симметрично положительных матриц с доминирующими диагональными
% элементами
k=0;
ok = false;
[m, n] = size(A);
x = zeros(m, 1);
if m == n && length(b) == m
    D = diag(diag(A));
    L = tril(A) - diag(diag(A));
    U = triu(A) - diag(diag(A));
    B = -D^(-1) * (L + U);
    g = D^(-1) * b;
    x_k = x0'; % previous 
    x_kplus1 = b'; % current
    q = normest(B); % rate of convergence %%скорость сходимости
    if normest(B)<1%(all(abs(eigs(B)) < 1))  % &&  all(diag(A) ~= 0))
        ok = true;
        while (k < kmax) && (max( abs(x_kplus1 - x_k) )/(1 - q) > epsilon)
            x_k = x_kplus1;
            for i = 1 : 1 : m
                sum = 0;
                for j = 1 : 1 : n
                    if (i ~= j)
                        sum = sum + A(i, j) * x_k(j);
                    end
                end
                x_kplus1(i) = (1/A(i, i)) * (b(i) - sum);
            end
            k = k + 1;
        end
    end
    x = x_kplus1';
end
end

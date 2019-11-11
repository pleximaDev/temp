function [x, ok, k] = my_successive_over_relaxation(A, b, x0, kmax, epsilon)
% итерационный вид
%(omega * L + D) * xk+1 = - (omega * U + (omega - 1) * D) * xk + omega * B
% для симметричных положительно определённых матриц
[m, n] = size(A);
D = diag(diag(A));
L = tril(A) - D;
U = triu(A) - D; 
P = (- (L + D) ^ (-1) ) * U;
k = 0;
x = x0;
T = eye(m) - D^(-1) * A;
p = max(abs(eigs(T))); % спектральный радиус матрицы
omega = 1 - (p/(1+sqrt(1-p^2)))^2;
ok = false;
if ((normest(P) < 1))
    ok = true;
    while ((k < kmax) && (norm(A * x - b) >= epsilon))
        
        for i = 1 : 1 : m
            Sum = 0;
            %{
            for j = 1 : 1 : n
                if i~=j
                    Sum = Sum + A(i, j) * x(j);
                end
                
            end
            %}
            for j = 1 : 1 : n
                if j < i
                    Sum = Sum + A(i, j) * x(j);
                end
                
            end
            for j = 1 : 1 : n
                if j > i
                    Sum = Sum + A(i, j) * x(j);
                end
                
            end
            x(i) = (1 - omega) * x(i) + (omega / A(i, i)) * (b(i) - Sum);
        end
        k = k + 1;
    end

end

end


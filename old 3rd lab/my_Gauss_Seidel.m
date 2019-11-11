function [x, ok, k] = my_Gauss_Seidel(A, b, x0, kmax, epsilon)
% итерационный вид
% (L + D) * xk+1 = -U * xk + B
% сходится быстрее, чем Якоби, если матрица симметричная и положительно
% определённая
[m, n] = size(A);
D = diag(diag(A));
L = tril(A) - D;

U = triu(A) - D;

P = - ( (L + D) ^ (-1) ) * U
k = 0;
x = x0;
C=zeros(m);
if ((normest(P) < 1))
    ok = true;
    for i = 1 : 1 : m
        for j = 1 : 1 : n
            if i~=j
                C(i,j)=-(A(i,j)/ A(i, i));
            end
        end
    end
    while ((k < kmax) && (normest(A * x - b) > epsilon))
        for i = 1 : 1 : n
            Sum = 0;
            d = b(i) / A(i, i);
            for j = 1:1:i-1
                    Sum = Sum + C(i, j) * x(j);
            end
            for j = i+1:1:n
                    Sum = Sum + C(i, j) * x(j);
            end
                x(i) = Sum + d;
        end
        k = k + 1;
    end
else
    ok = false;
    x = zeros(size(A, 1), 1);
end
end
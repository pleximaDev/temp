clc;
clear variables;
close all force;

a = 0;
b = 1;

Na = 8;
Nb = 10 * Na;

x1 = a : (b - a)/Na : b; %interpolated grid
x0 = a : (b - a)/Nb : b; %interpolation grid

N1 = length(x1);
N0 = length(x0);

v = 1; % frequency
omega = 2 * pi * v; 
f0 = interp_function(omega, x0); % analytic function model
f1 = interp_function(omega, x1); % experimental data model

[f2, x0] = nearest_neighbor(f1, x1, x0);
[f3, x0] = linear_interpolation(f1, x1, x0);
[f4, x0] = lagrange_polynomial(f1, x1, x0);
[f5, x0] = Newton_polynomial_w(f1, x1, x0, 'forward');
[f6, x0] = Newton_polynomial_w(f1, x1, x0, 'backward');
[f7, x0] = cubic_spline(f1, x1, x0);



clf;
subplot(2, 1, 1);

hold on;
grid on;
grid minor;

% plot(x0, f4);
plot(x0, f0, 'LineWidth', 1.5); 
stem(x1, f1, 'LineWidth', 1.5);
stairs(x0, f2, 'LineWidth', 1.5);
plot(x0, f3, 'LineWidth', 1.5);
plot(x0, f4, 'LineWidth', 1.5);

title('Interpolation');
legend('Experimental data model','Analytical model function','Nearest-neighbor interpolation','Linear interpolation','Lagrange');
legend('location','northeastoutside');
ylabel('f(x)');
xlabel('x');
hold off;

subplot(2, 1, 2);
hold on
plot(x0, f5, 'LineWidth', 1.5);
plot(x0, f6, 'LineWidth', 1.5);
plot(x0, f7, 'LineWidth', 1.5);
plot(x0, f0, 'LineWidth', 1.5);
stem(x1, f1, 'LineWidth', 1.5);
% hold on;
grid on;
grid minor;
% plot(x0, f6);
% plot(x0, f7);
% stem(x1, f1);
% plot(x0, f0);
legend('Newton polynomial forward','Newton polynomial backward','Cubic spline interpolation', 'Experimental data model','Analytical model function')
legend('location','northeastoutside');
ylabel('f(x)');
xlabel('x');
hold off;
% ylabel('f(x)');
% xlabel('x');
% hold off;




function f = interp_function(omega, x)
f = sin(omega * x);
end

function [f, x0] = nearest_neighbor(f1, x1, x0)
k = 1;
f = zeros(length(x0), 1);
for i = 1 : 1 : length(x0)
    if k < length(x1)
        delta1 = abs(x0(i) - x1(k));
        delta2 = abs(x0(i) - x1(k + 1));
        if delta2 <= delta1
            k = k + 1;
        end
    end
    f(i, 1) = f1(k);
end
end

function [f, x0] = linear_interpolation(f1, x1, x0)
k = 1;
f = zeros(length(x0), 1);
for i = 1 : 1 : length(x0)
    if x0(i) >= x1(k + 1) && k ~= length(x1) - 1
        k = k + 1;
    end
    f(i, 1) = f1(k) + ((f1(k + 1) - f1(k))/(x1(k + 1) - x1(k))) * (x0(i) - x1(k));
end

end

function [f, x0] = lagrange_polynomial(f1, x1, x0)
f = zeros(length(x0), 1);
shift = 1;
for l = 1 : 1 : length(x0)
    for k = 0 : 1 : length(x1) - shift
        basis_polynomial = 1;
        for i = 0 : 1 : length(x1) - shift
            if(i ~= k)
                basis_polynomial = basis_polynomial * (x0(l) - x1(i + shift))/(x1(k + shift) - x1(i + shift));
            end
        end
        f(l, 1) = f(l, 1) + basis_polynomial * f1(k + shift);
        
    end
end
end

%{
% function [f, x0] = Newton_polynomial(f1, x1, x0, Newton)
% f = zeros(length(x0), 1);
% finite_diff = 0;
% switch Newton
%     case 'forward'
%         for i = 1 : 1 : length(x0)
%             Polynomial = f1(1);
%             h = x1(1 + 1) - x1(1);
%             q = (x0(i) - f1(1))/h;
%             n = length(x1) - 1;
%             for j = 1 : 1 : n
%                 for v = 0 : 1 : j
%                     N = j;
%                     C = factorial(N)/(factorial(v) * factorial(N - v))
%                     finite_diff = finite_diff + (-1)^v * C * f1(1 + N - v)
%                 end
%                 Polynomial = Polynomial + (q - N + 1)/(factorial(N)) * finite_diff
%             end
%             
%             f(i) = Polynomial;
%         end
%         
% %         n = length(x1);
% %         for i=1 : 1 : length(x0);
% %             polynomial = f1(1);
% %             h = (x1(1 + 1) - x1(1));
% %             q=((x0(i)-x1(1))/h);
% %             for j=1:1:n
% % __________
% %             for v = 0 : 1 : j
% %                 N = j;
% %                 C = factorial(N)/(factorial(v) * factorial(N - v))
% %                 finite_diff = finite_diff + (-1)^v * C * f1(1 + N - v)
% %             end
% % __________
% %                 polynomial = polynomial + q(j)/factorial(j)*finite_diff;
% %                 q(j+1)=q(j)*(((x0(i)-x1(1))/h)-(j+1)+1);
% %             end
% %             f(i)=polynomial;
% %         end
%         
%     case 'backward'
%         for p=1:1:N0
%         end
% end
% end
%}



function [f, x0] = Newton_polynomial_w(f1, x1, x0, Newton)
f=zeros(length(x0), 1);
n = length(x1) - 1;

switch Newton
    case 'forward'
        for j = 1 : 1 : length(x0)
            polynomial=f1(1);
            h = (x1(2) - x1(1));
            q=((x0(j) - x1(1))/h);
            for i = 1 : 1 : n
                finite_diff = 0;
                for v = 0 : 1 : i
                    N = i;
                    C = factorial(N)/(factorial(v)*factorial(N - v));
                    finite_diff = finite_diff + (-1)^(v) * C *f1(1 + i - v);
                end
                polynomial = polynomial + q/factorial(i)*finite_diff;
                
                q = q * (((x0(j) - x1(1))/h) - (i + 1) + 1);
            end
            f(j) = polynomial;
        end
    case 'backward'
        for j = 1 : 1 : length(x0)
            polynomial=f1(n+1);  
            h = (x1(n+1)-x1(n));
            q=(x0(j)-x1(length(x1)))/h;
            for i = 1 : 1 : n
                finite_diff = 0;
                for v = 0 : 1 : i
                    N = i;
                    C = factorial(N)/(factorial(v)*factorial(N - v));
                    finite_diff = finite_diff + (-1)^(v)*C*f1((n + 1) - v);
                end
                polynomial = polynomial + q/factorial(i) * finite_diff;
                h = (x1(n+1)-x1(n));
                q = q * ((x0(j)-x1(length(x1)))/h + i + 1 - 1);
            end
            f(j)=polynomial;
        end
end

%     function y=divided_difference(n, Newton, f1, k)
%     y=0;
%     switch Newton
%         case 'forward'
%             for v=0:1:n
%                 y=y+(-1)^(v)*(factorial(n)/(factorial(v)*factorial((n-v))))*f1(k+n-v);
%             end
%         case 'backward'
%             for v=0:1:n
%                 y=y+(-1)^(v)*(factorial(n)/(factorial(v)*factorial((n-v))))*f1(k-v);
%             end
%     end
% 
%     end
end



% function [f, x0] = Newton_polynomial_test(f1, x1, x0, Newton)
% f = zeros(length(x0), 1);
% finite_diff = 0;
% switch Newton
%     case 'forward'
%         n = length(x1) - 1;
%         for p=1:1:N0
%             polNewton=f1(1);
%             h = (x1(2)-x1(1));
%             q=((x0(p)-x1(1))/h);
%             for i=1:1:n
%                 polNewton = polNewton + q/factorial(i)*divided_difference(i, Newton, f1, 1);
%                 h = (x1(2)-x1(1));
%                 q=q*(((x0(p)-x1(1))/h)-(i+1)+1);
%             end
%             f(p)=polNewton;
%         end
%         
%     case 'backward'
%         for p=1:1:N0
%         end
% end
% end

function [f, x0] = cubic_spline(f1, x1, x0)
shift = 1;
n = length(x1) - 1;
f = zeros(length(x0), 1);
c = 0;
c(n + 1) = 0;
K = zeros(length(x1), 1);
L = zeros(length(x1), 1);

for k = 2 : 1 : n
    h(k) =  x1(k + shift) - x1(k - 1 + shift);
    h(k - 1) =  x1(k - 1 + shift) - x1(k - 2 + shift);
    F = 3 * ((f1(k + 1) - f1(k))/h(k) - (f1(k)-f1(k-1))/h(k - 1));
    V=2 * (h(k) + h(k - 1));
    K(k) = (F - h(k - 1) * K(k - 1))/(V - h(k - 1) * L(k - 1));
    L(k) = h(k)/(V - h(k - 1) * L(k - 1));
end

for k = n : -1 : 2
    c(k) = K(k) - L(k) * c(k + 1);
end

for k = 1 : 1 : n
    h(k) = x1(k + shift) - x1(k - 1 + shift);
    d(k)=(c(k+1)-c(k))/(3*h(k));
    b(k)=(f1(k + shift) - f1(k - 1 + shift))/h(k)-c(k)*h(k)-d(k)*(h(k))^2;
    a(k)=f1(k - 1 + shift);
end

k = 1;
for i= 1 : 1 : length(x0)
    h(k) = x0(i) - x1(k);
    f(i)=a(k) + b(k) * h(k) + c(k) * (h(k))^2 + d(k) * (h(k))^3;
    
        if(x0(i) > x1(k + 1))
            k=k + 1;
        end
end 
% k = 2;
% for i = 1 : 1 : length(x0)
%     if x0(i) > x1(k) && k < n 
%         k = k + 1;
%     end
%     h = x0(i) - x1(k - 1);
% %     h(k) = x0(k) - x1(k - 1);
%     k
%     h
% %     s(k) = a(k) + b(k) * h(k) + c(k) * h(k)^2 + d(k) * h(k)^3;
%     s(i) = a(k) + b(k) * h + c(k) * h^2 + d(k) * h^3;
%     f(i) = s(i);
% end
end






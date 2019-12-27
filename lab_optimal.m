
clc;
clear variables;
close all force;

fnc_sca = @(x) (x + 1)^2 + 2;
fnc_vec = @(x) (x(1) + 1)^2 + (x(2) + 1)^2 + 2;

x = [4; 2];



% syms x y
% f = (x + 1)^2 + (y + 1)^2 + 2;
% hessian(f,[x,y])

a = -5;
b = 5;
eps = 1e-3;
Kmax = 1e6;

[Fmin_(1), Xmin_(1), k_(1)] = lab_optimal_sca(fnc_sca, a, b, eps, Kmax, "brute force");
[Fmin_(2), Xmin_(2), k_(2)] = lab_optimal_sca(fnc_sca, a, b, eps, Kmax, "Dichotomy");
[Fmin_(3), Xmin_(3), k_(3)] = lab_optimal_sca(fnc_sca, a, b, eps, Kmax, "Golden-section search");

methods = {'brute force'; 'Dichotomy';'Golden-section search'};
Fmin = Fmin_';
Xmin = Xmin_';
k = k_';

table(methods, Fmin, Xmin, k)

x0 = zeros(size(fnc_vec_fd(x, 0), 1), 1)

[Fmin2(1), Xmin2(:,1), k2(1)] = lab_optimal_vec(fnc_vec, @fnc_vec_fd, @fnc_vec_hessian, x0, eps, Kmax, "Pokoord");
[Fmin2(2), Xmin2(:,2), k2(2)] = lab_optimal_vec(fnc_vec, @fnc_vec_fd, @fnc_vec_hessian, x0, eps, Kmax, "Skoreish");
[Fmin2(3), Xmin2(:,3), k2(3)] = lab_optimal_vec(fnc_vec, @fnc_vec_fd, @fnc_vec_hessian, x0, eps, Kmax, "Newton");

methods = {'pokoord'; 'skoreish'; 'newton'};
Fmin = Fmin2';
Xmin = Xmin2';
k = k2';

table(methods, Fmin, Xmin, k)

%{
%\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\%
fmcOpt2 = optimoptions(@fmincon, 'Display','iter');
fmcOpt2.FiniteDifferenceType = 'central'
% fmcOpt2.OutputFcn = @hw_optimal_f;
% fmcOpt2.PlotFcns = @hw_optimal_f;
fmcOpt2.MaxFunctionEvaluations = 10000;
fmcOpt2.FiniteDifferenceStepSize = 1e-10;
fmcOpt2.TolX = 1e-9;
fmcOpt2.TolFun = 1e-9;
fmcOpt2.StepTolerance = 1e-9;
fmcOpt2.MaxIterations = 2000;
[xmin, fmin]=fmincon(fnc_vec,x0,[],[],[],[],[],[],[],fmcOpt2)
%\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\%
%}

%{
% fmcOpt = optimset('Display','iter-detailed','PlotFcns',@hw_optimal_f);
fmcOpt = optimset('Display','iter');
fmcOpt.MaxIter=2000;
fmcOpt.MaxFunEvals=10000;
fmcOpt.TolFun = 1e-9;
fmcOpt.TolX = 1e-9;
% fmcOpt.OutputFcn = @hw_optimal_f;
fmcOpt.FinDiffRelStep = 1e-10; %10^(-10)
fmcOpt.FinDiffType = 'central';
%_________________________________________________________________%
[xmin(:,1), fmin(1)] = fminsearch(fnc_sca, -5, fmcOpt)
%}


function [df] = fnc_vec_fd(x, i)
df1 = 2 * (x(1) + 1);
df2 = 2 * (x(2) + 1);
switch i 
    case 0
        df = [df1; df2];
    case 1
        df = df1;
    case 2
        df = df2;
    otherwise
        fprintf("You have chosen a nonexistent differential number i!");
end
end

function [H] = fnc_vec_hessian(x)
H = [2, 0; 0, 2];
end

function [Fmin, Xmin, k] = lab_optimal_sca(fnc, a, b, eps, Kmax, method)
k = 0;
switch method
    case "brute force"
        n = (b - a)/eps;
        xk = a;
        Xmin = xk;
        Fmin = fnc(Xmin);
        while(k < Kmax && xk < b)
            k = k + 1;
            xk = a + k * (b - a)/(n + 1);
            if(fnc(xk) < Fmin)
                Xmin = xk;
                Fmin = fnc(xk);
            end
%             Xmin = (fnc(xk) < Fmin) * xk + ~(fnc(xk) < Fmin) * Xmin;
%             Fmin = (fnc(xk) < Fmin) * fnc(xk) + ~(fnc(xk) < Fmin) * Fmin;
        end
        
    case 'Dichotomy'
        while(k < Kmax && abs(b - a) > eps)
            Xmin = (b + a)/2;
            k = k + 1;
            delta = (b - a)/4;
            if (fnc(Xmin - delta) >= fnc(Xmin + delta))
                a = Xmin;
            else
                b = Xmin;
            end
%             a = (fnc(Xmin - delta) >= fnc(Xmin + delta)) * Xmin + ~(fnc(Xmin - delta) >= fnc(Xmin + delta)) * a;
%             b = ~(fnc(Xmin - delta) >= fnc(Xmin + delta)) * Xmin + (fnc(Xmin - delta) >= fnc(Xmin + delta)) * b;
        end
        Xmin = (a + b)/2;
        Fmin = fnc(Xmin);
        
    case "Golden-section search"
%         phi = (1 + sqrt(5))/2;
        phi = 1.62;
        while(k < Kmax && abs(b - a) > eps)
            k = k + 1;
            Xmin = (b + a)/2;
            x1 = b - (b - a)/phi;
            x2 = a + (b - a)/phi;
            y1 = fnc(x1);
            y2 = fnc(x2);
            if(y1 >= y2)
                a = x1;
                x1 = x2;
                x2 = a + (b - a)/phi;
            else
                b = x2;
                x2 = x1;
                x1 = b - (b - a)/phi;
            end
        end
        Xmin = (a + b)/2;
        Fmin = fnc(Xmin);
    otherwise
        fprintf("You have chosen a nonexistent method!");
        
end

end



function [Fmin, Xmin, k] = lab_optimal_vec(fnc, fncd, fncH, x0, eps, Kmax, method)
k = 0;
switch method 
    case "Pokoord"
        % xkn = xk
        %xk = xk_1
        xk = x0;
        xk_1 = [5; 5];
        while abs(fnc(xk) - fnc(xk_1)) >= eps && norm(xk - xk_1) >= eps && k < Kmax
            k = k + 1;
            xk_1 = xk;
            for i = 1 : 1 : length(xk)
                xk(i) = xk_1(i) + lambda(xk_1, eps, Kmax, i) * (-fncd(xk_1,i));
            end
            Xmin = xk;
            Fmin = fnc(Xmin);
        end
        
    case "Skoreish"
        xk = x0;
        xk_1 = [5; 5];
        while abs(fnc(xk) - fnc(xk_1)) >= eps && norm(xk - xk_1) >= eps && k < Kmax
            xk_1 = xk;
            k = k + 1;
            xk = xk_1 + lambda(xk_1, eps, Kmax, 0) * (-fncd(xk_1, 0));
        end
        Xmin = xk;
        Fmin = fnc(Xmin);
        
    case "Newton"
        xk = x0;
        xk_1 = [5; 5];
        while abs(fnc(xk) - fnc(xk_1)) >= eps && norm(xk - xk_1) >= eps && k < Kmax
            xk_1 = xk;
            k = k + 1;
            xk = xk_1 + (fncH(xk)^(-1)) * (-fncd(xk_1, 0));
        end
        Xmin = xk;
        Fmin = fnc(Xmin);
    otherwise
        fprintf("You have chosen a nonexistent method!");
end

function [y] = lambda(xk, eps, Kmax, i)
a = 0;
b = 2;
j = 0;
%         phi = (1 + sqrt(5))/2;
phi = 1.62;
while(j < Kmax && abs(b - a) > eps)
    j = j + 1;
    Xmin = (b + a)/2;
    x1 = b - (b - a)/phi;
    x2 = a + (b - a)/phi;
    y1 = fnc(xk + x1 * (-fncd(xk, i)) );
    y2 = fnc(xk + x2 * (-fncd(xk, i)) );
    if(y1 >= y2)
        a = x1;
        x1 = x2;
        x2 = a + (b - a)/phi;
    else
        b = x2;
        x2 = x1;
        x1 = b - (b - a)/phi;
    end
end
Xmin = (a + b)/2;
% Fmin = fnc(Xmin);
% Fmin = fnc(xk + Xmin * fncd(xk, i));

y = Xmin;
end

end




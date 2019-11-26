clc;
clear variables;
close all force;
format short;

fnc = @(t) lab_diff_f(t);

% d - derivative order
% p - finite diff order
% [A, C, b, divider] = C_coeff(d, p, method)
% divider of whole finite diff == 1/divider * (factorial(d)/1)

% rats(num2str(C(1)))

% [A, C, b, divider, d, p] = C_coeff(1, 4, "backward")
% [str] = str_finite_diff(C, d, p, divider, "backward")

a = 0.2;
b = 0.7;
n = 20;



[finite_diff_df(:, 1, 1)] = lab_diff_do(fnc, a, b, n, 1, "forward"); 
[finite_diff_df(:, 2, 1)] = lab_diff_do(fnc, a, b, n, 1, "backward");
[finite_diff_df(:, 1, 2)] = lab_diff_do(fnc, a, b, n, 2, "forward");
[finite_diff_df(:, 2, 2)] = lab_diff_do(fnc, a, b, n, 2, "backward");
[finite_diff_df(:, 3, 2)] = lab_diff_do(fnc, a, b, n, 2, "central");
[finite_diff_df(:, 1, 3)] = lab_diff_do(fnc, a, b, n, 4, "forward");
[finite_diff_df(:, 2, 3)] = lab_diff_do(fnc, a, b, n, 4, "backward");
[finite_diff_df(:, 3, 3)] = lab_diff_do(fnc, a, b, n, 4, "central");
[finite_diff_df(:, 1, 4)] = lab_diff_do(fnc, a, b, n, 6, "forward");
[finite_diff_df(:, 2, 4)] = lab_diff_do(fnc, a, b, n, 6, "backward");
[finite_diff_df(:, 3, 4)] = lab_diff_do(fnc, a, b, n, 6, "central");




% #2
%-------------
h = (b - a)/n;
t = a : h : b;
%-------------
% 
[f] = lab_diff_f(t);
[df] = lab_diff_df(t);
% 

%============================plot============================
figure('Name', 'Finite differences for n = 20','Numbertitle', 'off')
clf
subplot(2, 2, 1);
plot(t, finite_diff_df(:, 1, 1), 'LineWidth', 1.5);
title("First order finite differences, except central")
grid on;
grid minor;

hold on;
plot(t, finite_diff_df(:, 2, 1), 'LineWidth', 1.5);
hold on;
plot(t, df, 'LineWidth', 1.5);
legend('Forward finite difference', 'Backward finite difference', ...
    'Default derivative');
hold off;

subplot(2, 2, 2);
plot(t, finite_diff_df(:, 1, 2), 'LineWidth', 1.5);
title("Second order finite differences")
grid on;
grid minor;
hold on;
plot(t, finite_diff_df(:, 2, 2), 'LineWidth', 1.5);
hold on;
plot(t, df, 'LineWidth', 1.5);
hold on;
plot(t, finite_diff_df(:, 3, 2), 'LineWidth', 1.5);
legend('Forward finite difference','Backward finite difference', ...
  'Default derivative', 'Central finite difference');
hold off;

subplot(2, 2, 3);
plot(t, finite_diff_df(:, 1, 3), 'LineWidth', 1.5);
title("Fourth order finite differences")
grid on;
grid minor;
hold on;
plot(t, finite_diff_df(:, 2, 3), 'LineWidth', 1.5);
hold on;
plot(t, df, 'LineWidth', 1.5);
hold on;
plot(t, finite_diff_df(:, 3, 3), 'LineWidth', 1.5);
legend('Forward finite difference','Backward finite difference', ...
  'Default derivative', 'Central finite difference');
hold off;

subplot(2, 2, 4);
plot(t, finite_diff_df(:, 1, 4), 'LineWidth', 1.5);
title("Sixth order finite differences")
grid on;
grid minor;
hold on;
plot(t, finite_diff_df(:, 2, 4), 'LineWidth', 1.5);
hold on;
plot(t, df, 'LineWidth', 1.5);
hold on;
plot(t, finite_diff_df(:, 3, 4), 'LineWidth', 1.5);
legend('Forward finite difference','Backward finite difference', ...
  'Default derivative', 'Central finite difference');
hold off;
%============================================================


% 1 order
error(1, 1) = mean(abs(~isnan(finite_diff_df(:, 1, 1)) - df')./abs(df)');
error(2, 1) = mean(abs(~isnan(finite_diff_df(:, 2, 1)) - df')./abs(df)');
error(3, 1) = 0;

% 2 order
error(1, 2) = mean(abs(~isnan(finite_diff_df(:, 1, 2)) - df')./abs(df)');
error(2, 2) = mean(abs(~isnan(finite_diff_df(:, 2, 2)) - df')./abs(df)');
error(3, 2) = mean(abs(~isnan(finite_diff_df(:, 3, 2)) - df')./abs(df)');

% 4 order
error(1, 3) = mean(abs(~isnan(finite_diff_df(:, 1, 3)) - df')./abs(df)');
error(2, 3) = mean(abs(~isnan(finite_diff_df(:, 2, 3)) - df')./abs(df)');
error(3, 3) = mean(abs(~isnan(finite_diff_df(:, 3, 3)) - df')./abs(df)');

% 6 order
error(1, 4) = mean(abs(~isnan(finite_diff_df(:, 1, 4)) - df')./abs(df)');
error(2, 4) = mean(abs(~isnan(finite_diff_df(:, 2, 4)) - df')./abs(df)');
error(3, 4) = mean(abs(~isnan(finite_diff_df(:, 3, 4)) - df')./abs(df)');

%==========================bar==========================
figure('Name', 'Error for n = 20','Numbertitle', 'off')
clf
bar(error);
title('Error estimation');
legend('First order','Second order','Fourth order','Sixth order', 'Location', 'NorthEastOutside');
ax = gca;
ax.XTickLabel = {'Forward','Backward','Central'};
grid on;
grid minor;
%=======================================================


% 2.5


n = 1000;
%-------------
h = (b - a)/n;
t = a : h : b;
%-------------

[df] = lab_diff_df(t);
[finite_diff_df_1000(:, 1, 1)] = lab_diff_do(fnc, a, b, n, 1, "forward"); 
[finite_diff_df_1000(:, 2, 1)] = lab_diff_do(fnc, a, b, n, 1, "backward");
[finite_diff_df_1000(:, 1, 2)] = lab_diff_do(fnc, a, b, n, 2, "forward");
[finite_diff_df_1000(:, 2, 2)] = lab_diff_do(fnc, a, b, n, 2, "backward");
[finite_diff_df_1000(:, 3, 2)] = lab_diff_do(fnc, a, b, n, 2, "central");
[finite_diff_df_1000(:, 1, 3)] = lab_diff_do(fnc, a, b, n, 4, "forward");
[finite_diff_df_1000(:, 2, 3)] = lab_diff_do(fnc, a, b, n, 4, "backward");
[finite_diff_df_1000(:, 3, 3)] = lab_diff_do(fnc, a, b, n, 4, "central");
[finite_diff_df_1000(:, 1, 4)] = lab_diff_do(fnc, a, b, n, 6, "forward");
[finite_diff_df_1000(:, 2, 4)] = lab_diff_do(fnc, a, b, n, 6, "backward");
[finite_diff_df_1000(:, 3, 4)] = lab_diff_do(fnc, a, b, n, 6, "central");



%============================plot============================
figure('Name', 'Finite differences for n = 1000','Numbertitle', 'off')
clf
subplot(2, 2, 1);
plot(t, finite_diff_df_1000(:, 1, 1), 'LineWidth', 1.5);
title("First order finite differences, except central")
grid on;
grid minor;

hold on;
plot(t, finite_diff_df_1000(:, 2, 1), 'LineWidth', 1.5);
hold on;
plot(t, df, 'LineWidth', 1.5);
legend('Forward finite difference', 'Backward finite difference', ...
    'Default derivative');
hold off;

subplot(2, 2, 2);
plot(t, finite_diff_df_1000(:, 1, 2), 'LineWidth', 1.5);
title("Second order finite differences")
grid on;
grid minor;
hold on;
plot(t, finite_diff_df_1000(:, 2, 2), 'LineWidth', 1.5);
hold on;
plot(t, df, 'LineWidth', 1.5);
hold on;
plot(t, finite_diff_df_1000(:, 3, 2), 'LineWidth', 1.5);
legend('Forward finite difference','Backward finite difference', ...
  'Default derivative', 'Central finite difference');
hold off;

subplot(2, 2, 3);
plot(t, finite_diff_df_1000(:, 1, 3), 'LineWidth', 1.5);
title("Fourth order finite differences")
grid on;
grid minor;
hold on;
plot(t, finite_diff_df_1000(:, 2, 3), 'LineWidth', 1.5);
hold on;
plot(t, df, 'LineWidth', 1.5);
hold on;
plot(t, finite_diff_df_1000(:, 3, 3), 'LineWidth', 1.5);
legend('Forward finite difference','Backward finite difference', ...
  'Default derivative', 'Central finite difference');
hold off;

subplot(2, 2, 4);
plot(t, finite_diff_df_1000(:, 1, 4), 'LineWidth', 1.5);
title("Sixth order finite differences")
grid on;
grid minor;
hold on;
plot(t, finite_diff_df_1000(:, 2, 4), 'LineWidth', 1.5);
hold on;
plot(t, df, 'LineWidth', 1.5);
hold on;
plot(t, finite_diff_df_1000(:, 3, 4), 'LineWidth', 1.5);
legend('Forward finite difference','Backward finite difference', ...
  'Default derivative', 'Central finite difference');
hold off;
%============================================================





%==========================bar==========================
figure('Name', 'Error for n = 1000','Numbertitle', 'off')
clf
bar(error);
title('Error estimation');
legend('First order','Second order','Fourth order','Sixth order', 'Location', 'NorthEastOutside');
ax = gca;
ax.XTickLabel = {'Forward','Backward','Central'};
grid on;
grid minor;
%=======================================================

syms x1  x2
J = jacobian([(x1 + 1)^2 + (x2 + 1)^2 + 2, (x2 - 1)^2 + (x2 - x1)^2], [x1, x2])
% j = horzcat([f1; f2], J)



% syms x1 x2
f1(x1, x2) = (x1 + 1)^2 + (x2 + 1)^2 + 2;
f2(x1, x2) = (x2 - 1)^2 + (x2 - x1)^2;
f = [f1, f2];
f(1, 1)
(6 + 1)^2 + (2 + 1)^2 + 2
(6 - 1)^2 + (6 - 2)^2


[    2*x1 + 2,        2*x2 + 2]
[ 2*x1 - 2*x2, 4*x2 - 2*x1 - 2]

df1_dx1 = 2 * x1 + 2;
df2_dx1 = 2 * x1 - 2 * x2;
df1_dx2 = 2 * x2 + 2;
df2_dx2 = 4 * x2 - 2 * x1 - 2;
my_Jacobian(x1, x2) = [df1_dx1, df1_dx2; df2_dx1, df2_dx2]
my_Jacobian(1,1)
% f = f(1, 1)
% subs(f)



% 
% function [J] = Jacobi(x1, x2)
% end


function [f] = lab_diff_f(t)
v = 4; % [Hz]
omega = 2 * pi * v; % [rad/s]
f = 1.16 * t + 0.13 * sin(omega * t) - 0.89 * t.^2;
end

function [df] = lab_diff_df(t)
v = 4; % [Hz]
omega = 2 * pi * v; % [rad/s]
df = 1.16 + 0.13 * omega * cos(omega * t) - 1.78 * t;
end

function [df, t] = lab_diff_do(fnc, a, b, n, k, method)
h = (b - a)/n;
t = a : h : b;
df(1 : 1 : length(t)) = NaN;
switch k 
    case 1
        switch method
            
                case 'forward'
                    for i = 1 : 1 : length(t) - 1
                        df(i) = (fnc(t(i) + h) - fnc(t(i)))/h;
                    end
                case 'backward'
                    for i = 1 + 1 : 1 : length(t)
                        df(i) = (fnc(t(i)) - fnc(t(i) - h))/h;
                    end
                case 'central'
                    % Nope
                otherwise
                    fprintf("Error occured while entering method's name.");
        end
    case 2
        switch method
            case 'forward'
                for i = 1 : 1 : length(t) - 2
                    df(i) = ( -fnc(t(i) + 2*h) + 4 * fnc(t(i) + h) - 3 * fnc(t(i)))/(2 * h);
                end
            case 'backward'
                for i = 1 + 2 : 1 : length(t)
                    df(i) = ( 3 * fnc(t(i)) - 4 * fnc(t(i) - h) + fnc(t(i) - 2*h))/(2 * h);
                end
            case 'central'
                for i = 1 + 1  : 1 : length(t) - 1
                    df(i) = (fnc(t(i) + h) - fnc(t(i) - h))/(2 * h);
                end
            otherwise
                fprintf("Error occured while entering method's name.");
        end
    case 4
        switch method
            case 'forward'
                for i = 1 : 1 : length(t) - 4
                    df(i) = ((-3) * fnc(t(i) + 4*h) + 16 * fnc(t(i) + 3*h) - 36 * fnc(t(i) + 2*h) + 48 * fnc(t(i) + h) - 25 * fnc(t(i)))/(12 * h);
                end
            case 'backward'
                for i = 1 + 4 : 1 : length(t)
                    df(i) = (25 * fnc(t(i)) - 48 * fnc(t(i) - h) + 36 * fnc(t(i) - 2 * h) - 16 * fnc(t(i) - 3 * h) + 3 * fnc(t(i) - 4 * h))/(12*h);
                end
            case 'central'
                for i = 1 + 2 : 1 : length(t) - 2
                    df(i) = (-fnc(t(i) + 2*h) + 8 * fnc(t(i) + h) - 8 * fnc(t(i) - h) + fnc(t(i) - 2 * h))/(12 * h);
                end
            otherwise
                fprintf("Error occured while entering method's name.");
        end
    case 6
        switch method
            case 'forward'
                for i = 1 : 1 : length(t) - 6
                    df(i) = ((-1/6) * fnc(t(i) + 6*h) + (6/5) * fnc(t(i) + 5*h) - (15/4) * fnc(t(i) + 4*h) + (20/3) * fnc(t(i) + 3*h) - (15/2) * fnc(t(i) + 2*h) + 6 * fnc(t(i) + 1 * h) -(49/20)* fnc(t(i)))/h;
                end
            case 'backward'
                for i = 1 + 6 : 1 : length(t)
                  df(i) = ((49/20) * fnc(t(i)) - 6 * fnc(t(i) - h) + (15/2) * fnc(t(i) - 2 * h) + (-20/3) * fnc(t(i) - 3 * h) + (15/4) * fnc(t(i) - 4 * h) - (6/5) * fnc(t(i) - 5 * h) + (1/6)* fnc(t(i) - 6 * h))/h;
                end
            case 'central'
                for i = 1 + 3 : 1 : length(t) - 3
                    df(i) = ((-1/60) * fnc(t(i) - 3*h) + (3/20)*fnc(t(i) - 2*h) - (3/4) * fnc(t(i) - h) + ...
                     (3/4)* fnc(t(i) + h) + (-3/20)* fnc(t(i) + 2*h) + (1/60) * fnc(t(i) + 3*h))/h;
                end
            otherwise
                fprintf("Error occured while entering method's name.");
        end
end
end
















function [A, C, b, divider, d, p] = C_coeff(d, p, method)
format rat;
switch method
    case "forward"
        i_min = 0;
        i_max = d + p - 1;
        n = size(i_min : 1 : i_max, 2);
        i = (0 : 1 : n);
        for t = 1 : 1 : n
            A(t, :) =  [(i_min : i_max).^i(t)];
        end
        b = zeros(n, 1);
        b(d + 1) = 1;
        C = linsolve(A, b);
        [N, D] = rat(C);
        divider = max(D);
        C = C * divider;
%         A = [(i_min : i_max).^0; (i_min : i_max).^1; (i_min : i_max).^2; (i_min : i_max).^3]

    case "backward"
        i_min = -(d + p - 1);
        i_max = 0;
        n = size(i_min : 1 : i_max, 2);
        i = (0 : 1 : n);
        for t = 1 : 1 : n
            A(t, :) =  [(i_min : i_max).^i(t)];
        end
        b = zeros(n, 1);
        b(d + 1) = 1;
        C = linsolve(A, b);
        [N, D] = rat(C);
        divider = max(D);
        C = C * divider;
    case "centered"
        i_min = -fix((d + p - 1)/2); % i_max = -i_min
        i_max = fix((d + p - 1)/2);
        n = size(i_min : 1 : i_max, 2);
        i = (0 : 1 : n);
        for t = 1 : 1 : n
            A(t, :) =  [(i_min : 1 : i_max).^i(t)];
        end
        b = zeros(n, 1);
        b(d + 1) = 1;
        C = linsolve(A, b);
        C = linsolve(A, b);
        [N, D] = rat(C);
        divider = max(D);
        C = C * divider;
    otherwise
        fprintf("You have chosen a nonexistent method!");
end
% format short;
end

function [str] = str_finite_diff(C, d, p, divider, method)
switch method
    case "forward"
        
%         C = C(end:-1:1)
        if (d ~= 1)
            start = ['$$F^' int2str(d) '(x) = '];
        else
            start = ['$$F(x) = '];
        end
        if (p >= 6)
            format rat;
            C
            C = C/divider
            C_str = rats(C);
            format short;
            if(C(end) ~= 1)
                str = ['{ ' C_str(end, :) ' * F(x + ' int2str(size(C, 1) - 1) ' * h) + '];
    %             str = ['{ ' int2str(C(end)) ' * F(x + ' int2str(p) ' * h) + '];
            else
                str = ['{F(x + ' int2str(size(C, 1) - 1) ' * h) + '];
            end
            for i = size(C, 1) - 2 : -1 : 1 %i = p - 1 : -1 : 1
                if(C(i+1) < -1e-11)
                    str(end-1:end) = [];
                end
                if(C(i+1) > 1e-11 || C(i+1) < -1e-11)
                    if(abs(C(i+1)) ~= 1)
                        str = [str C_str(i+1, :) ' * F(x + ' int2str(i) ' * h) + '];
    %                     str = [str int2str(C(i+1)) ' * F(x + ' int2str(i) ' * h) + '];
                    elseif(C(i+1) == 1)
                        str = [str ' F(x + ' int2str(i) ' * h) + '];
                    else
                        str = [str ' - F(x + ' int2str(i) ' * h) + '];
                    end
                end
            end
            if(C(1) < -1e-11)
                    str(end-1:end) = [];
            end
            if(abs(C(1)) ~= 1)
                str = [str C_str(1, :) ' * F(x)'];
    %             str = [str int2str(C(1)) ' * F(x)'];
            elseif(C(1) == 1)
                str = [str ' F(x)'];
            else
                str = [str ' - F(x)'];
            end
        else
            %________________________________p<6________________________________%
            if(C(end) ~= 1)
%                 str = ['{ ' C_str(end, :) ' * F(x + ' int2str(p) ' * h) + '];
                str = ['{ ' int2str(C(end)) ' * F(x + ' int2str(size(C, 1)-1) ' * h) + '];
            else
                str = ['{F(x + ' int2str(size(C, 1) - 1) ' * h) + '];
            end
            for i = size(C, 1) - 2 : -1 : 1
                if(C(i+1) < -1e-11)
                    str(end-1:end) = [];
                end
                if(C(i+1) > 1e-11 || C(i+1) < -1e-11)
                    if(abs(C(i+1)) ~= 1)
%                         str = [str C_str(i+1, :) ' * F(x + ' int2str(i) ' * h) + '];
                        str = [str int2str(C(i+1)) ' * F(x + ' int2str(i) ' * h) + '];
                    elseif(C(i+1) == 1)
                        str = [str ' F(x + ' int2str(i) ' * h) + '];
                    else
                        str = [str ' - F(x + ' int2str(i) ' * h) + '];
                    end
                end
            end
            if(C(1) < -1e-11)
                    str(end-1:end) = [];
            end
            if(abs(C(1)) ~= 1)
%                 str = [str C_str(1) ' * F(x)'];
                str = [str int2str(C(1)) ' * F(x)'];
            elseif(C(1) == 1)
                str = [str ' F(x)'];
            else
                str = [str ' - F(x)'];
            end
        end
        
        if(d >= 1)
            if(d == 1)
                if(p < 6)
                    str = [str '\over' int2str(divider) ' * h}'];
                else
                    str = [str '\over h}'];
                end
            else
                if(p < 6)
                    str = [str '\over' int2str(factorial(d)/divider) ' * h^' int2str(d) '}'];
                else
                    str = [str '\over h^' int2str(d) '}'];
                end
            end
        end
        str = [str '$$'];
        start = [start str];
        str = start;
        %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    case "backward"
        C = C(end:-1:1)
        if (d ~= 1)
            start = ['$$F^' int2str(d) '(x) = '];
        else
            start = ['$$F(x) = '];
        end
        
        if (p >= 6)
            format rat;
            C
            C = C/divider
            C_str = rats(C);
            format short;
        
            if(abs(C(1)) ~= 1)
                str = ['{ ' C_str(1, :) ' * F(x) + '];
%                 str = ['{ ' int2str(C(1)) ' * F(x) + '];
            else
                str = ['{F(x) +'];
            end

            if(C(2) > 0)
                str = [str C_str(C(2, :)) ' * F(x - h) + '];
%                 str = [str int2str(C(2)) ' * F(x - h) + '];
            elseif (C(2) < 0)
                str(end-1:end) = [];
                str = [str C_str(2, :) ' * F(x - h) + '];
%                 str = [str int2str(C(2)) ' * F(x - h) + '];
            end
            for i = 3 : 1 : size(C, 1)
                if(C(i) < -1e-11)
                    str(end-1:end) = [];
                end
                if(C(i) > 1e-11 || C(i) < -1e-11)
                    if(abs(C(i)) ~= 1)
                        str = [str C_str(i, :) ' * F(x - ' int2str(i) - 1 ' * h) + '];
%                         str = [str int2str(C(i)) ' * F(x - ' int2str(i) - 1 ' * h) + '];
                    elseif(C(i) == 1)
                        str = [str ' F(x - ' int2str(i) - 1 ' * h) + '];
                    else
                        str = [str ' - F(x - ' int2str(i) - 1 ' * h) + '];
                    end
                end
            end
           
            %____p<6______
        else
            if(abs(C(1)) ~= 1)
                str = ['{ ' int2str(C(1)) ' * F(x) + '];
            else
                str = ['{F(x) +'];
            end

            if(C(2) > 0)
                str = [str int2str(C(2)) ' * F(x - h) + '];
            elseif (C(2) < 0)
                str(end-1:end) = [];
                str = [str int2str(C(2)) ' * F(x - h) + '];
            end
            for i = 3 : 1 : size(C, 1)
                if(C(i) < -1e-11)
                    str(end-1:end) = [];
                end
                if(C(i) > 1e-11 || C(i) < -1e-11)
                    if(abs(C(i)) ~= 1)
                        str = [str int2str(C(i)) ' * F(x - ' int2str(i) - 1 ' * h) + '];
                    elseif(C(i) == 1)
                        str = [str ' F(x - ' int2str(i) - 1 ' * h) + '];
                    else
                        str = [str ' - F(x - ' int2str(i) - 1 ' * h) + '];
                    end
                end
            end
            
        end
        %_____________________________
            
        str(end-1:end) = [];
        if(d >= 1)
            if(d == 1)
                if(p < 6)
                    str = [str '\over' int2str(divider) ' * h}'];
                else
                    str = [str '\over h}'];
                end
            else
                if(p < 6)
                    str = [str '\over' int2str(factorial(d)/divider) ' * h^' int2str(d) '}'];
                else
                    str = [str '\over h^' int2str(d) '}'];
                end
            end
        end
        str = [str '$$'];
        start = [start str];
        str = start;
    case "centered"
        if (d ~= 1)
            start = ['$$F^' int2str(d) '(x) = '];
        else
            start = ['$$F(x) = '];
        end
        str = '{';
        i = 1;
        if(p >= 6)
            format rat;
            C = C/divider
            C_str = rats(C)
            format short;
            for j = -fix(size(C, 1)/2) : 1 : fix(size(C, 1)/2)
                if(C(i) > 1e-11 || C(i) < -1e-11)
                    if (j < 0)
                        if(abs(C(i)) ~= 1)
                            if(C(i) > 1e-11)
                                str = [str C_str(i, :) ' * F(x ' int2str(j) ' * h) + '];
                            elseif((C(i) < -1e-11) && (i > 1))

                                str(end-1:end) = [];
                                str = [str C_str(i, :) ' * F(x ' int2str(j) ' * h) + '];
                            else
                                str = [str C_str(i, :) ' * F(x ' int2str(j) ' * h) + '];
                            end

        %                         str = [str int2str(C(i)) ' * F(x - ' int2str(i) - 1 ' * h) + '];
                        elseif(C(i) == 1)
                            str = [str ' F(x ' int2str(j) ' * h) + '];
                        else
                            str = [str ' - F(x ' int2str(j) ' * h) + '];
                        end
                    elseif(j > 0)
                        if(abs(C(i)) ~= 1)
                            if(C(i) > 1e-11)
                                str = [str C_str(i, :) ' * F(x + ' int2str(j) ' * h) + '];
                            elseif((C(i) < -1e-11) && (i > 1))
                                str(end-1:end) = [];
                                str = [str C_str(i, :) ' * F(x + ' int2str(j) ' * h) + '];
                            else
                                str = [str C_str(i, :) ' * F(x + ' int2str(j) ' * h) + '];
                            end

        %                         str = [str int2str(C(i)) ' * F(x - ' int2str(i) - 1 ' * h) + '];
                        elseif(C(i) == 1)
                            str = [str ' F(x + ' int2str(j) ' * h) + '];
                        else
                            str(end-1:end) = [];
                            str = [str ' - F(x + ' int2str(j) ' * h) + '];
                        end
                    else
                        if(abs(C(i)) ~= 1)
                            str = [str C_str(i, :) ' * F(x) + '];

        %                         str = [str int2str(C(i)) ' * F(x - ' int2str(i) - 1 ' * h) + '];
                        elseif(C(i) == 1)
                            str = [str ' F(x) + '];
                        else
                            str = [str ' - F(x) + '];
                        end
                    end

                end
                i = i + 1;
            end
        else %(p < 6)
            for j = -fix(size(C, 1)/2) : 1 : fix(size(C, 1)/2)
                if(C(i) > 1e-11 || C(i) < -1e-11)
                    if (j < 0)
                        if(abs(C(i)) ~= 1)
                            if(C(i) > 1e-11)
                                str = [str int2str(C(i)) ' * F(x ' int2str(j) ' * h) + '];
                            elseif((C(i) < -1e-11) && (i > 1))

                                str(end-1:end) = [];
                                str = [str int2str(C(i)) ' * F(x ' int2str(j) ' * h) + '];
                            else
                                str = [str int2str(C(i)) ' * F(x ' int2str(j) ' * h) + '];
                            end

        %                         str = [str int2str(C(i)) ' * F(x - ' int2str(i) - 1 ' * h) + '];
                        elseif(C(i) == 1)
                            str = [str ' F(x ' int2str(j) ' * h) + '];
                        else
                            str = [str ' - F(x ' int2str(j) ' * h) + '];
                        end
                    elseif(j > 0)
                        if(abs(C(i)) ~= 1)
                            if(C(i) > 1e-11)
                                str = [str int2str(C(i)) ' * F(x + ' int2str(j) ' * h) + '];
                            elseif((C(i) < -1e-11) && (i > 1))
                                str(end-1:end) = [];
                                str = [str int2str(C(i)) ' * F(x + ' int2str(j) ' * h) + '];
                            else
                                str = [str int2str(C(i)) ' * F(x + ' int2str(j) ' * h) + '];
                            end

        %                         str = [str int2str(C(i)) ' * F(x - ' int2str(i) - 1 ' * h) + '];
                        elseif(C(i) == 1)
                            str = [str ' F(x + ' int2str(j) ' * h) + '];
                        else
                            str(end-1:end) = [];
                            str = [str ' - F(x + ' int2str(j) ' * h) + '];
                        end
                    else
                        if(abs(C(i)) ~= 1)
                            str = [str int2str(C(i)) ' * F(x) + '];

        %                         str = [str int2str(C(i)) ' * F(x - ' int2str(i) - 1 ' * h) + '];
                        elseif(C(i) == 1)
                            str = [str ' F(x) + '];
                        else
                            str = [str ' - F(x) + '];
                        end
                    end

                end
                i = i + 1;
            end
        end
        %++++++++++++++++++++++++++++++++
        str(end-1:end) = [];
        if(d >= 1)
            if(d == 1)
                if(p < 6)
                    str = [str '\over' int2str(divider) ' * h}'];
                else
                    str = [str '\over h}'];
                end
            else
                if(p < 6)
                    str = [str '\over' int2str(factorial(d)/divider) ' * h^' int2str(d) '}'];
                else
                    str = [str '\over h^' int2str(d) '}'];
                end
            end
        end
        str = [str '$$'];
        start = [start str];
        str = start;
end
figure('Name', 'Approximation of derivative', 'NumberTitle', 'off', ...
    'Position', [200 300 1400 300]); % [left bottom width height]
clf
plot(1, 1);
title([int2str(d) ' order derivative approximated by ' int2str(p) ...
    ' order ' char(method) ' finite difference']);
ylim([0, 3]);
xlim([0, 8]);
text(0.5, 1.5, str, 'Interpreter', 'Latex', 'Fontsize', 16);
end

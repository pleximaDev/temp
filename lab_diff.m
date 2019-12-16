clc;
clear variables;
close all force;
format short;

% d - derivative order
% p - finite diff order
% [A, C, b, divider] = C_coeff(d, p, method)
% divider of whole finite diff == 1/divider * (factorial(d)/1)

[A, C, b, divider, d, p] = C_coeff(1, 4, "backward");
[str] = str_finite_diff(C, d, p, divider, "backward");





for n = 20 : 1000 - 20 : 1000
    fnc = @(t) lab_diff_f(t);
    a = 0.2;
    b = 0.7;
    n

    %-------------
    h = (b - a)/n;
    t = a : h : b;
    %-------------

    [f] = lab_diff_f(t);
    [df] = lab_diff_df(t);
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




    %============================plot============================
    Name = ['Finite differences for n = ', char(string(n))];
    figure('Name', Name,'Numbertitle', 'off')

    clf
    subplot(2, 2, 1);
    plot(t, finite_diff_df(:, 1, 1));
    title("First order finite differences, except central")
    grid on;
    grid minor;

    hold on;
    plot(t, finite_diff_df(:, 2, 1));
    hold on;
    plot(t, df);
    legend('Forward finite difference', 'Backward finite difference', ...
        'Default derivative');
    hold off;

    subplot(2, 2, 2);
    plot(t, finite_diff_df(:, 1, 2) )
    title("Second order finite differences")
    grid on;
    grid minor;
    hold on;
    plot(t, finite_diff_df(:, 2, 2) )
    hold on;
    plot(t, df);
    hold on;
    plot(t, finite_diff_df(:, 3, 2) )
    legend('Forward finite difference','Backward finite difference', ...
      'Default derivative', 'Central finite difference');
    hold off;

    subplot(2, 2, 3);
    plot(t, finite_diff_df(:, 1, 3))
    title("Fourth order finite differences")
    grid on;
    grid minor;
    hold on;
    plot(t, finite_diff_df(:, 2, 3))
    hold on;
    plot(t, df);
    hold on;
    plot(t, finite_diff_df(:, 3, 3))
    legend('Forward finite difference','Backward finite difference', ...
      'Default derivative', 'Central finite difference');
    hold off;

    subplot(2, 2, 4);
    plot(t, finite_diff_df(:, 1, 4))
    title("Sixth order finite differences")
    grid on;
    grid minor;
    hold on;
    plot(t, finite_diff_df(:, 2, 4))
    hold on;
    plot(t, df);
    hold on;
    plot(t, finite_diff_df(:, 3, 4))
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
    Name = ['Error for n = ', char(string(n))];
    figure('Name', Name,'Numbertitle', 'off')
    clf
    bar(error);
    title('Error estimation');
    legend('First order','Second order','Fourth order','Sixth order', 'Location', 'NorthEastOutside');
    ax = gca;
    ax.XTickLabel = {'Forward','Backward','Central'};
    grid on;
    grid minor;
    %=======================================================
    clear variables;
end



% % % % #2
% % % %-------------
% % % h = (b - a)/n;
% % % t = a : h : b;
% % % %-------------
% % % % 
% % % [f] = lab_diff_f(t);
% % % [df] = lab_diff_df(t);
% % % % 
% % % 
% % % %============================plot============================
% % % figure('Name', 'Finite differences for n = 20','Numbertitle', 'off')
% % % clf
% % % subplot(2, 2, 1);
% % % plot(t, finite_diff_df(:, 1, 1));
% % % title("First order finite differences, except central")
% % % grid on;
% % % grid minor;
% % % 
% % % hold on;
% % % plot(t, finite_diff_df(:, 2, 1));
% % % hold on;
% % % plot(t, df);
% % % legend('Forward finite difference', 'Backward finite difference', ...
% % %     'Default derivative');
% % % hold off;
% % % 
% % % subplot(2, 2, 2);
% % % plot(t, finite_diff_df(:, 1, 2) )
% % % title("Second order finite differences")
% % % grid on;
% % % grid minor;
% % % hold on;
% % % plot(t, finite_diff_df(:, 2, 2) )
% % % hold on;
% % % plot(t, df);
% % % hold on;
% % % plot(t, finite_diff_df(:, 3, 2) )
% % % legend('Forward finite difference','Backward finite difference', ...
% % %   'Default derivative', 'Central finite difference');
% % % hold off;
% % % 
% % % subplot(2, 2, 3);
% % % plot(t, finite_diff_df(:, 1, 3))
% % % title("Fourth order finite differences")
% % % grid on;
% % % grid minor;
% % % hold on;
% % % plot(t, finite_diff_df(:, 2, 3))
% % % hold on;
% % % plot(t, df);
% % % hold on;
% % % plot(t, finite_diff_df(:, 3, 3))
% % % legend('Forward finite difference','Backward finite difference', ...
% % %   'Default derivative', 'Central finite difference');
% % % hold off;
% % % 
% % % subplot(2, 2, 4);
% % % plot(t, finite_diff_df(:, 1, 4))
% % % title("Sixth order finite differences")
% % % grid on;
% % % grid minor;
% % % hold on;
% % % plot(t, finite_diff_df(:, 2, 4))
% % % hold on;
% % % plot(t, df);
% % % hold on;
% % % plot(t, finite_diff_df(:, 3, 4))
% % % legend('Forward finite difference','Backward finite difference', ...
% % %   'Default derivative', 'Central finite difference');
% % % hold off;
% % % %============================================================
% % % 
% % % % 1 order
% % % error(1, 1) = mean(abs(~isnan(finite_diff_df(:, 1, 1)) - df')./abs(df)');
% % % error(2, 1) = mean(abs(~isnan(finite_diff_df(:, 2, 1)) - df')./abs(df)');
% % % error(3, 1) = 0;
% % % 
% % % % 2 order
% % % error(1, 2) = mean(abs(~isnan(finite_diff_df(:, 1, 2)) - df')./abs(df)');
% % % error(2, 2) = mean(abs(~isnan(finite_diff_df(:, 2, 2)) - df')./abs(df)');
% % % error(3, 2) = mean(abs(~isnan(finite_diff_df(:, 3, 2)) - df')./abs(df)');
% % % 
% % % % 4 order
% % % error(1, 3) = mean(abs(~isnan(finite_diff_df(:, 1, 3)) - df')./abs(df)');
% % % error(2, 3) = mean(abs(~isnan(finite_diff_df(:, 2, 3)) - df')./abs(df)');
% % % error(3, 3) = mean(abs(~isnan(finite_diff_df(:, 3, 3)) - df')./abs(df)');
% % % 
% % % % 6 order
% % % error(1, 4) = mean(abs(~isnan(finite_diff_df(:, 1, 4)) - df')./abs(df)');
% % % error(2, 4) = mean(abs(~isnan(finite_diff_df(:, 2, 4)) - df')./abs(df)');
% % % error(3, 4) = mean(abs(~isnan(finite_diff_df(:, 3, 4)) - df')./abs(df)');
% % % 
% % % %==========================bar==========================
% % % figure('Name', 'Error for n = 20','Numbertitle', 'off')
% % % clf
% % % bar(error);
% % % title('Error estimation');
% % % legend('First order','Second order','Fourth order','Sixth order', 'Location', 'NorthEastOutside');
% % % ax = gca;
% % % ax.XTickLabel = {'Forward','Backward','Central'};
% % % grid on;
% % % grid minor;
% % % %=======================================================


% % % % % % 2.5
% % % % % 
% % % % % 
% % % % % n = 1000;
% % % % % %-------------
% % % % % h = (b - a)/n;
% % % % % t = a : h : b;
% % % % % %-------------
% % % % % 
% % % % % [df] = lab_diff_df(t);
% % % % % [finite_diff_df_1000(:, 1, 1)] = lab_diff_do(fnc, a, b, n, 1, "forward"); 
% % % % % [finite_diff_df_1000(:, 2, 1)] = lab_diff_do(fnc, a, b, n, 1, "backward");
% % % % % [finite_diff_df_1000(:, 1, 2)] = lab_diff_do(fnc, a, b, n, 2, "forward");
% % % % % [finite_diff_df_1000(:, 2, 2)] = lab_diff_do(fnc, a, b, n, 2, "backward");
% % % % % [finite_diff_df_1000(:, 3, 2)] = lab_diff_do(fnc, a, b, n, 2, "central");
% % % % % [finite_diff_df_1000(:, 1, 3)] = lab_diff_do(fnc, a, b, n, 4, "forward");
% % % % % [finite_diff_df_1000(:, 2, 3)] = lab_diff_do(fnc, a, b, n, 4, "backward");
% % % % % [finite_diff_df_1000(:, 3, 3)] = lab_diff_do(fnc, a, b, n, 4, "central");
% % % % % [finite_diff_df_1000(:, 1, 4)] = lab_diff_do(fnc, a, b, n, 6, "forward");
% % % % % [finite_diff_df_1000(:, 2, 4)] = lab_diff_do(fnc, a, b, n, 6, "backward");
% % % % % [finite_diff_df_1000(:, 3, 4)] = lab_diff_do(fnc, a, b, n, 6, "central");
% % % % % 
% % % % % 
% % % % % 
% % % % % %============================plot============================
% % % % % figure('Name', 'Finite differences for n = 1000','Numbertitle', 'off')
% % % % % clf
% % % % % subplot(2, 2, 1);
% % % % % plot(t, finite_diff_df_1000(:, 1, 1));
% % % % % title("First order finite differences, except central")
% % % % % grid on;
% % % % % grid minor;
% % % % % 
% % % % % hold on;
% % % % % plot(t, finite_diff_df_1000(:, 2, 1));
% % % % % hold on;
% % % % % plot(t, df);
% % % % % legend('Forward finite difference', 'Backward finite difference', ...
% % % % %     'Default derivative');
% % % % % hold off;
% % % % % 
% % % % % subplot(2, 2, 2);
% % % % % plot(t, finite_diff_df_1000(:, 1, 2) )
% % % % % title("Second order finite differences")
% % % % % grid on;
% % % % % grid minor;
% % % % % hold on;
% % % % % plot(t, finite_diff_df_1000(:, 2, 2) )
% % % % % hold on;
% % % % % plot(t, df);
% % % % % hold on;
% % % % % plot(t, finite_diff_df_1000(:, 3, 2) )
% % % % % legend('Forward finite difference','Backward finite difference', ...
% % % % %   'Default derivative', 'Central finite difference');
% % % % % hold off;
% % % % % 
% % % % % subplot(2, 2, 3);
% % % % % plot(t, finite_diff_df_1000(:, 1, 3))
% % % % % title("Fourth order finite differences")
% % % % % grid on;
% % % % % grid minor;
% % % % % hold on;
% % % % % plot(t, finite_diff_df_1000(:, 2, 3))
% % % % % hold on;
% % % % % plot(t, df);
% % % % % hold on;
% % % % % plot(t, finite_diff_df_1000(:, 3, 3))
% % % % % legend('Forward finite difference','Backward finite difference', ...
% % % % %   'Default derivative', 'Central finite difference');
% % % % % hold off;
% % % % % 
% % % % % subplot(2, 2, 4);
% % % % % plot(t, finite_diff_df_1000(:, 1, 4))
% % % % % title("Sixth order finite differences")
% % % % % grid on;
% % % % % grid minor;
% % % % % hold on;
% % % % % plot(t, finite_diff_df_1000(:, 2, 4))
% % % % % hold on;
% % % % % plot(t, df);
% % % % % hold on;
% % % % % plot(t, finite_diff_df_1000(:, 3, 4))
% % % % % legend('Forward finite difference','Backward finite difference', ...
% % % % %   'Default derivative', 'Central finite difference');
% % % % % hold off;
% % % % % %============================================================
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % % 1 order
% % % % % error(1, 1) = mean(abs(~isnan(finite_diff_df_1000(:, 1, 1)) - df')./abs(df)');
% % % % % error(2, 1) = mean(abs(~isnan(finite_diff_df_1000(:, 2, 1)) - df')./abs(df)');
% % % % % error(3, 1) = 0;
% % % % % 
% % % % % % 2 order
% % % % % error(1, 2) = mean(abs(~isnan(finite_diff_df_1000(:, 1, 2)) - df')./abs(df)');
% % % % % error(2, 2) = mean(abs(~isnan(finite_diff_df_1000(:, 2, 2)) - df')./abs(df)');
% % % % % error(3, 2) = mean(abs(~isnan(finite_diff_df_1000(:, 3, 2)) - df')./abs(df)');
% % % % % 
% % % % % % 4 order
% % % % % error(1, 3) = mean(abs(~isnan(finite_diff_df_1000(:, 1, 3)) - df')./abs(df)');
% % % % % error(2, 3) = mean(abs(~isnan(finite_diff_df_1000(:, 2, 3)) - df')./abs(df)');
% % % % % error(3, 3) = mean(abs(~isnan(finite_diff_df_1000(:, 3, 3)) - df')./abs(df)');
% % % % % 
% % % % % % 6 order
% % % % % error(1, 4) = mean(abs(~isnan(finite_diff_df_1000(:, 1, 4)) - df')./abs(df)');
% % % % % error(2, 4) = mean(abs(~isnan(finite_diff_df_1000(:, 2, 4)) - df')./abs(df)');
% % % % % error(3, 4) = mean(abs(~isnan(finite_diff_df_1000(:, 3, 4)) - df')./abs(df)');
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % %==========================bar==========================
% % % % % figure('Name', 'Error for n = 1000','Numbertitle', 'off')
% % % % % clf
% % % % % bar(error);
% % % % % title('Error estimation');
% % % % % legend('First order','Second order','Fourth order','Sixth order', 'Location', 'NorthEastOutside');
% % % % % ax = gca;
% % % % % ax.XTickLabel = {'Forward','Backward','Central'};
% % % % % grid on;
% % % % % grid minor;
% % % % % %=======================================================




syms x1  x2
% J = jacobian([(x1 + 1)^2 + (x2 + 1)^2 + 2, (x2 - 1)^2 + (x2 - x1)^2], [x1, x2])

% % % f1(x1, x2) = (x1 + 1)^2 + (x2 + 1)^2 + 2;
% % % f2(x1, x2) = (x2 - 1)^2 + (x2 - x1)^2;
% % % f = [f1, f2];

f1 = @(x1,x2) (x1 + 1)^2 + (x2 + 1)^2 + 2;
f2 = @(x1,x2) (x2 - 1)^2 + (x2 - x1)^2;
% % % f1 = @(x) (x(1) + 1)^2 + (x(2) + 1)^2 + 2;
% % % f2 = @(x) (x(2) - 1)^2 + (x(2) - x(1))^2;
f = {f1,f2};


% [    2*x1 + 2,        2*x2 + 2]
% [ 2*x1 - 2*x2, 4*x2 - 2*x1 - 2]

% % % % 2 * x1 + 2          = 1
% % % % 2 * x2 + 2          = 0
% % % % 2 * x1 - 2 * x2     = 0
% % % % 4 * x2 - 2 * x1 - 2 = 1




close all;
x = [1,1]
df1_dx1 = 2 * x1 + 2;
df2_dx1 = 2 * x1 - 2 * x2;
df1_dx2 = 2 * x2 + 2;
df2_dx2 = 4 * x2 - 2 * x1 - 2;

my_Jacobian(x1, x2) = [df1_dx1, df1_dx2; df2_dx1, df2_dx2];
my_Jacobian_1 = {@(x1,x2) 2 * x1 + 2; @(x1,x2) 2 * x1 - 2 * x2; ...
    @(x1,x2) 2 * x2 + 2; @(x1,x2) 4 * x2 - 2 * x1 - 2}


df1_dx1 =@(x1,x2) 2 * x1 + 2;
df2_dx1 =@(x1,x2) 2 * x1 - 2 * x2;
df1_dx2 =@(x1,x2) 2 * x2 + 2;
df2_dx2 =@(x1,x2) 4 * x2 - 2 * x1 - 2;
my_Jacobian_2= {df1_dx1, df1_dx2; df2_dx1, df2_dx2};



x = [1,1]
[J_forward] = Jacobi_finite_diff(x, f, "forward")
[J_backward] = Jacobi_finite_diff(x, f, "backward")
[J_central] = Jacobi_finite_diff(x, f, "central")
my = my_Jacobian(x(1), x(2))

% A = abs((J1 - my_Jacobian(x(1), x(2)))./my_Jacobian(x(1), x(2)))
analytical_Jacobian = double(my_Jacobian(x(1), x(2)));
numerator = (J_forward - analytical_Jacobian);
divider = analytical_Jacobian;
error_Jacobi(1) = norm(abs(numerator), 'fro')/norm(abs(divider), 'fro');

numerator = (J_backward - analytical_Jacobian);
divider = analytical_Jacobian;
error_Jacobi(2) = norm(abs(numerator), 'fro')/norm(abs(divider), 'fro');

numerator = (J_central - analytical_Jacobian);
divider = analytical_Jacobian;
error_Jacobi(3) = norm(abs(numerator), 'fro')/norm(abs(divider), 'fro');


figure('Name', 'Relative error of Jacobian calculation','Numbertitle', 'off')
clf
bar(error_Jacobi)
title('Error estimation');
ax = gca;
ax.XTickLabel = {'Forward','Backward','Central'};
grid on;
grid minor;

clc

my_Jacobian(1,1)

format short
[J] = Broyden_funct(x, @f12, 3)

return


f
% [J] = Broyden([1,1], f, 3)
[J_0] = Jacobi_finite_diff([1,2], f, "forward")

[J] = Broyden_cell(J_0, x, f, 3)

% f1(x1, x2) = (x1 + 1)^2 + (x2 + 1)^2 + 2;
% f2(x1, x2) = (x2 - 1)^2 + (x2 - x1)^2;
% f = [f1, f2]
% f(1, 1)



function [f3] = f12(x)
f1 = (x(1) + 1)^2 + (x(2) + 1)^2 + 2;
f2 = (x(2) - 1)^2 + (x(2) - x(1))^2;
f3 = [f1,f2];
end

function [J] = Jacobi_finite_diff(x, f, method)
syms x1 x2
J = NaN;
h = 0.01;
switch method
    case 'forward'
        for i = 1 : 1 : length(x)
            for j = 1 : 1 : length(x)
                J(j, i) = (f{j}(x(1) + h * (2 - i), x(2) + h * (i - 1)) - ...
                    f{j}(x(1), x(2)))/h;
            end
        end
    case 'backward'
        for i = 1 : 1 : length(x)
            for j = 1 : 1 : length(x)
                J(j, i) = (f{j}(x(1), x(2)) - ...
                    f{j}(x(1) - h * (2 - i), x(2) - h * (i - 1)))/h;
            end
        end
    case 'central'
        for i = 1 : 1 : length(x)
            for j = 1 : 1 : length(x)
                J(j, i) = (f{j}(x(1) + h * (2 - i), x(2) + h * (i - 1)) - ...
                    f{j}(x(1) - h * (2 - i), x(2) - h * (i - 1)))/(2 * h);
            end
        end
    otherwise
        fprintf("You have chosen a nonexistent method!");
end
end



function [J] = Jacobi_finite_diff_syms(x, f, method)
syms x1 x2
J = NaN;
h = 0.01;
fbody = formula(f); % (x(1), x(2)))
f1(x1,x2) = fbody(1);
f2(x1,x2) = fbody(2);
switch method
    case 'forward'
        for i = 1 : 1 : length(x)
            for j = 1 : 1 : length(x)
                J(j, i) = (subs(fbody(j), {x1,x2}, {x(1) + h * (2 - i), x(2) + h * (i - 1)}) - ...
                    subs(fbody(j), {x1,x2}, {x(1), x(2)}))/h;
            end
        end
    case 'backward'
        for i = 1 : 1 : length(x)
            for j = 1 : 1 : length(x)
                J(j, i) = (subs(fbody(j), {x1,x2}, {x(1), x(2)}) - ...
                    subs(fbody(j), {x1,x2}, {x(1) - h * (2 - i), x(2) - h * (i - 1)}))/h;
            end
        end
    case 'central'
        for i = 1 : 1 : length(x)
            for j = 1 : 1 : length(x)
                J(j, i) = (subs(fbody(j), {x1,x2}, {x(1) + h * (2 - i), x(2) + h * (i - 1)}) - ...
                    subs(fbody(j), {x1,x2}, {x(1) - h * (2 - i), x(2) - h * (i - 1)}))/(2 * h);
            end
        end
    otherwise
        fprintf("You have chosen a nonexistent method!");
end
end


function [J] = Broyden_funct(x, f12, kmax)
x_k_1 = [1, 1];
x_k = [2, 2];
J = eye(2);
for i = 1 : 1 : kmax
%         x_k = x_k_1 - (J^(-1) * f12(x_k_1)')';
        J = J + ((f12(x_k) - f12(x_k_1) ...
        - (J * (x_k - x_k_1)')')./((x_k - x_k_1) * (x_k - x_k_1)')) * (x_k - x_k_1)'
        x_k_1 = x_k;
end
end


function [J] = Broyden_cell(J_0, x, f, kmax)
% % % % % J = eye(length(x));
% J = [ 6,  4; 2, -2];
% % % % x_prev = [2,1]
% % % % % % % % % % % % % % % % x_prev = x;
% % % % % % % % % % % % % % % % J = J_0;
x_prev = [1, 1]
J = ones(2)


% f_x = [f{1}(x(1),x(2)), f{2}(x(1),x(2))]
f_prev = [f{1}(x_prev(1),x_prev(2)), f{2}(x_prev(1),x_prev(2))]
% x_cur = x_prev -  (J^(-1) * f_prev')'

for i = 1 : 1 : kmax
    for j = 1 : 1 : length(x)
        f_prev = [f{1}(x_prev(1),x_prev(2)), f{2}(x_prev(1),x_prev(2))]
%         x_cur = x_prev - (J^(-1) * f_prev')';
        x_cur = x_prev + 0.01;
        J 
        J^(-1)
        f_cur = [f{1}(x_cur(1),x_cur(2)), f{2}(x_cur(1),x_cur(2))]
        ((f_cur - f_prev - (J * (x_cur - x_prev)')')./((x_cur - x_prev) * (x_cur - x_prev)')) * (x_cur - x_prev)'
        J = J + ((f_cur - f_prev ...
        - (J * (x_cur - x_prev)')')./((x_cur - x_prev) * (x_cur - x_prev)')) * (x_cur - x_prev)'
%     return
        x_prev = x_cur;
    end
end
end

function [J] = Broyden(x, f, kmax)
J = 2 * eye(length(x));
h = 1e-2;
x_prev = x;
x_cur = x_prev - f(x(1),x(2)) * J^(-1)
for i = 1 : 1 : kmax
    for j = 1 : 1 : length(x)
        J = J + ((f(x_cur(1), x_cur(2)) - f(x_prev(1), x_prev(2)) ...
        - J * (x_cur - x_prev)')/(norm(x_cur - x_prev))^2) * (x_cur - x_prev)';
        x_prev = x_cur;
        x_cur = x_prev - f(x_prev(1), x_prev(2)) * J^(-1);
    end
end
end


%{
function [J] = Jacobi_finite_diff_by_element(x, f, method)
syms x1 x2
J = NaN;
h = 0.01;
fbody = formula(f); % (x(1), x(2)))
f1(x1,x2) = fbody(1);
f2(x1,x2) = fbody(2);
switch method
    case 'forward'
        for i = 1 : 1 : length(x)
            J(1, i) = (f1(x(1) + h * (2 - i), x(2) + h * (i - 1)) - f1(x(1), x(2)))/h;
            J(2, i) = (f2(x(1) + h * (2 - i), x(2) + h * (i - 1)) - f2(x(1), x(2)))/h;
        end
        
%             J(1, 1) = (f1(x(1) + h, x(2)) - f1(x(1), x(2)))/h;
%             J(2, 1) = (f2(x(1) + h, x(2)) - f2(x(1), x(2)))/h;
%             J(1, 2) = (f1(x(1), x(2) + h) - f1(x(1), x(2)))/h;
%             J(2, 2) = (f2(x(1), x(2) + h) - f2(x(1), x(2)))/h;
    case 'backward'
    case 'central'
    otherwise
        fprintf("You have chosen a nonexistent method!");
end
end
%}


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

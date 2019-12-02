clc;
clear variables;
close all force;
format short;

% d - derivative order
% p - finite diff order
% [A, C, b, divider] = C_coeff(d, p, method)
% divider of whole finite diff == 1/divider * (factorial(d)/1)

% [A, C, b, divider, d, p] = C_coeff(1, 4, "backward")
% [str] = str_finite_diff(C, d, p, divider, "backward")

fnc = @(t) lab_diff_f(t);

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
plot(t, finite_diff_df_1000(:, 1, 1));
title("First order finite differences, except central")
grid on;
grid minor;

hold on;
plot(t, finite_diff_df_1000(:, 2, 1));
hold on;
plot(t, df);
legend('Forward finite difference', 'Backward finite difference', ...
    'Default derivative');
hold off;

subplot(2, 2, 2);
plot(t, finite_diff_df_1000(:, 1, 2) )
title("Second order finite differences")
grid on;
grid minor;
hold on;
plot(t, finite_diff_df_1000(:, 2, 2) )
hold on;
plot(t, df);
hold on;
plot(t, finite_diff_df_1000(:, 3, 2) )
legend('Forward finite difference','Backward finite difference', ...
  'Default derivative', 'Central finite difference');
hold off;

subplot(2, 2, 3);
plot(t, finite_diff_df_1000(:, 1, 3))
title("Fourth order finite differences")
grid on;
grid minor;
hold on;
plot(t, finite_diff_df_1000(:, 2, 3))
hold on;
plot(t, df);
hold on;
plot(t, finite_diff_df_1000(:, 3, 3))
legend('Forward finite difference','Backward finite difference', ...
  'Default derivative', 'Central finite difference');
hold off;

subplot(2, 2, 4);
plot(t, finite_diff_df_1000(:, 1, 4))
title("Sixth order finite differences")
grid on;
grid minor;
hold on;
plot(t, finite_diff_df_1000(:, 2, 4))
hold on;
plot(t, df);
hold on;
plot(t, finite_diff_df_1000(:, 3, 4))
legend('Forward finite difference','Backward finite difference', ...
  'Default derivative', 'Central finite difference');
hold off;
%============================================================












% 1 order
error(1, 1) = mean(abs(~isnan(finite_diff_df_1000(:, 1, 1)) - df')./abs(df)');
error(2, 1) = mean(abs(~isnan(finite_diff_df_1000(:, 2, 1)) - df')./abs(df)');
error(3, 1) = 0;

% 2 order
error(1, 2) = mean(abs(~isnan(finite_diff_df_1000(:, 1, 2)) - df')./abs(df)');
error(2, 2) = mean(abs(~isnan(finite_diff_df_1000(:, 2, 2)) - df')./abs(df)');
error(3, 2) = mean(abs(~isnan(finite_diff_df_1000(:, 3, 2)) - df')./abs(df)');

% 4 order
error(1, 3) = mean(abs(~isnan(finite_diff_df_1000(:, 1, 3)) - df')./abs(df)');
error(2, 3) = mean(abs(~isnan(finite_diff_df_1000(:, 2, 3)) - df')./abs(df)');
error(3, 3) = mean(abs(~isnan(finite_diff_df_1000(:, 3, 3)) - df')./abs(df)');

% 6 order
error(1, 4) = mean(abs(~isnan(finite_diff_df_1000(:, 1, 4)) - df')./abs(df)');
error(2, 4) = mean(abs(~isnan(finite_diff_df_1000(:, 2, 4)) - df')./abs(df)');
error(3, 4) = mean(abs(~isnan(finite_diff_df_1000(:, 3, 4)) - df')./abs(df)');




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
f(1, 3)


% [    2*x1 + 2,        2*x2 + 2]
% [ 2*x1 - 2*x2, 4*x2 - 2*x1 - 2]

close all;
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

% return

x = [1,1]
fbody = formula(f(x(1), x(2)))
f1 = fbody(1)
f2 = fbody(2)
f1(1,1)

fbody2 = formula(f(x1,x2))
lel = fbody2(1)
lel(1,1)
whos fbody2
subs(fbody2(1)/3, {x1,x2}, {1, 2})
% return

clc
format short
x = [1; 1]
[J] = Jacobi_finite_diff(x, f, "forward")
my_Jacobian(x(1), x(2))
J = jacobian([(x1 + 1)^2 + (x2 + 1)^2 + 2, (x2 - 1)^2 + (x2 - x1)^2], [x1,x2])

subs(J(1),{1,1})


function [J] = Jacobi_finite_diff(x, f, method)
syms x1 x2
J = NaN;
h = 0.01;
fbody = formula(f); % (x(1), x(2)))
f1(x1,x2) = fbody(1);
f2(x1,x2) = fbody(2);

% df(i) = (fnc(t(i) + h) - fnc(t(i)))/h;
switch method
    case 'forward'
        for i = 1 : 1 : length(x)
            for j = 1 : 1 : length(x)
                J(j, i) = (subs(fbody(j), {x1,x2}, {x(1) + h * (2 - i), x(2) + h * (i - 1)}) - ...
                    subs(fbody(j), {x1,x2}, {x(1), x(2)}))/h
            end
        end
    case 'backward'
    case 'central'
    otherwise
        fprintf("You have chosen a nonexistent method!");
end
end

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




function [f] = lab_diff_f(t)
v = 4; % [Hz]
omega = 2 * pi * v; % [rad/s]
f = 1.16 * t + 0.13 * sin(omega * t) - 0.89 * t.^2;
end

function [df] = lab_diff_df(t)
v = 4 % [Hz]
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



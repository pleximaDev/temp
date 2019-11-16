clc;
clear variables;
close all force;


fnc = @(t) lab_diff_f(t);

% d - derivative order
% p - finite diff order
% [A, C, b, divider] = C_coeff(d, p, method)
% divider of whole finite diff == 1/divider * (factorial(d)/1)

[A, C, b, divider, d, p] = C_coeff(3, 4, "forward")

[A, C, b, divider, d, p] = C_coeff(3, 4, "centered")
% C = C(end:-1:1)
% 1/8
clc
[A, C, b, divider, d, p] = C_coeff(1, 6, "forward")


C_str = rats(C)

% rats(num2str(C(1)))
clc
% [A, C, b, divider, d, p] = C_coeff(1, 8, "backward")
% [str] = str_finite_diff(C, d, p, divider, "backward")

a = 0.2
b = 0.7
n = 20

[df, t] = lab_diff_do(fnc, a, b, n, 1, "forward");
[df, t] = lab_diff_do(fnc, a, b, n, 1, "backward");
[df, t] = lab_diff_do(fnc, a, b, n, 2, "forward");
[df, t] = lab_diff_do(fnc, a, b, n, 2, "backward");
[df, t] = lab_diff_do(fnc, a, b, n, 2, "central");
[df, t] = lab_diff_do(fnc, a, b, n, 4, "forward");
[df, t] = lab_diff_do(fnc, a, b, n, 4, "backward");
[df, t] = lab_diff_do(fnc, a, b, n, 4, "central");
[df, t] = lab_diff_do(fnc, a, b, n, 6, "forward");
[df, t] = lab_diff_do(fnc, a, b, n, 6, "backward");
[df, t] = lab_diff_do(fnc, a, b, n, 6, "central");
df = df'

% #2
%-------------
h = (b - a)/n;
t = a : h : b;
%-------------

[f] = lab_diff_f(t);



[df2] = lab_diff_df(t);
df2 = df2'



figure(1)
clf
plot(t, df)
hold on
plot(t, df2)
plot(t, f)
grid on 
grid minor
hold off



function [f] = lab_diff_f(t)
v = 4 % [Hz]
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
switch k 
    case 1
        switch method
            case "forward"
                df = (fnc(t + h) - fnc(t))/h;
            case "backward"
                df = (fnc(t) - fnc(t - h))/h;
            case "central"
                % Nope
            otherwise
                fprintf("Error occured while entering method's name.");
        end
    case 2
        switch method
            case "forward"
                df = ( -fnc(t + 2*h) + 4 * fnc(t + h) - 3 * fnc(t))/(2 * h);
            case "backward"
                df = ( 3 * fnc(t) - 4 * fnc(t - h) + fnc(t - 2*h))/(2 * h);
            case "central"
                df = (fnc(t + h) - fnc(t - h))/(2 * h);
            otherwise
                fprintf("Error occured while entering method's name.");
        end
    case 4
        switch method
            case "forward"
                df = ((-3) * fnc(t + 4*h) + 16 * fnc(t + 3*h) - 36 * fnc(t + 2*h) + 48 * fnc(t + h) - 25 * fnc(t))/(12 * h);
            case "backward"
                df = (25 * fnc(t) - 48 * fnc(t - h) + 36 * fnc(t - 2 * h) - 16 * fnc(t - 3 * h) + 3 * fnc(t - 4 * h))/(12*h);
            case "central"
                df = (-fnc(t + 2*h) + 8 * fnc(t + h) - 8 * fnc(t - h) + fnc(t - 2 * h))/(12 * h);
            otherwise
                fprintf("Error occured while entering method's name.");
        end
    case 6
        switch method
            case "forward"
                df = ((-1/6) * fnc(t + 6*h) + (6/5) * fnc(t + 5*h) - (15/4) * fnc(t + 4*h) + (20/3) * fnc(t + 3*h) - (15/2) * fnc(t + 2*h) + 6 * fnc(t + 1 * h) -(49/20)* fnc(t))/h;
            case "backward"
                df = ((49/20) * fnc(t) - 6 * fnc(t - h) + (15/2) * fnc(t - 2 * h) + (-20/3) * fnc(t - 3 * h) + (15/4) * fnc(t - 4 * h) - (6/5) * fnc(t - 5 * h) + (1/6)* fnc(t - 6 * h))/h;
            case "central"
                df = ((-1/60) * fnc(t - 3*h) + (3/20)*fnc(t - 2*h) - (3/4) * fnc(t - h) + ...
                    (3/4)* fnc(t + h) + (-3/20)* fnc(t + 2*h) + (1/60) * fnc(t + 3*h))/h;
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


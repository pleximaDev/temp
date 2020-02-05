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

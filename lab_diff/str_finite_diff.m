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

clc;
clear variables;
close all force;


% график итераций неправильный, k должны быть целыми, применимость Якоби


A = [1 -0.2589 -0.3093; -0.2589 1 -0.2705; -0.3093 -0.2705 1]
b = [1; 1; 1]
x0 = zeros(3, 1)
epsilon = 1e-10
kmax = 1000



2.2

[x, ok, k] = my_Jacobi(A, b, x0, kmax, epsilon)


[x, ok, k] = my_Gauss_Seidel(A, b, x0, kmax, epsilon)


[x, ok, k] = my_successive_over_relaxation(A, b, x0, kmax, epsilon)






K = 16;                             
b=randn(K, 1);
D=cell(4,1);

% 1) Matrica A poloj. opred. s dominir diag elem
A0 = randn(K);
I = eye(K, K);
A0 = tril(A0);

A = A0 * A0' + 5 * K * I;
D{1}={A,b};

% 2) matrica A simmietrichnaya, otricatelno opredelennaya, s domin diag elem
A0 = randn(K);
A0 = tril(A0);

A = A0 * A0' - 5 * K * I;
D{2}={A,b};

% 3) matrica A Neumann - razrejennaya otric opred matreye

A = gallery('neumann', K);
D{3}={A,b};
% 4) matrica A Neumamn razrejennaya otric opred matr v matlabe v vide POLNOY MATRICI(full)

A = full(A);
D{4}={A,b};

save('lab_slau_data.mat', 'D', '-v7');

clear variables;

load('lab_slau_data.mat');


% 3.1 
sizeA=16;
x0=zeros(sizeA,1);
X=zeros(sizeA, 3); 
T=zeros(4, 3);
K=zeros(1,3);
kmax=1000;
epsilon = 1e-10;

for i = 1 : 1 : 4
A = D{i}{1};
b = D{i}{2};
tic
[x, ok, k]=my_Jacobi(A, b, x0, kmax, epsilon);
t=toc;

T(i, 1)=t*ok;
X(:, 1, i) = x;

tic
[x, ok, k] = my_Gauss_Seidel(A, b, x0, kmax, epsilon);
t=toc;

T(i, 2)=t*ok;
X(:, 2, i) = x;

tic
[x, ok, k] = my_successive_over_relaxation(A, b, x0, kmax, epsilon);
t=toc;

T(i, 3)=t*ok;
X(:, 3, i) = x;
end


N = 1000;
T=zeros(4, 3);

for i = 1 : 1 : 4
%%%%%%%%%%%%%%%%%%%%% Jacobi %%%%%%%%%%%%%%%%%%%%%
    A = D{i}{1};
    b = D{i}{2};
    timeVector=zeros(N, 1);
    
    for j = 1 : 1 : N
        tic
    [x, ok, k]=my_Jacobi(A, b, x0, kmax, epsilon);
    timeVector(j, 1)=toc*ok;
    
    if ~ok 
        break
    end   

    end
    T(i, 1)= mean(timeVector); %zapicivaem srednee znachenie v T
    K(i, 1)=k*ok;


  
%%%%%%%%%%%%%%%%%%%%% Gauss-Seidel %%%%%%%%%%%%%%%%%%%%%


    A = D{i}{1};
    b = D{i}{2};
    timeVector=zeros(N, 1);
    
    for j = 1 : 1 : N
        tic
    [x, ok, k] = my_Gauss_Seidel(A, b, x0, kmax, epsilon);
    timeVector(j, 1)=toc*ok;
    if ~ok 
        break
    end   
     %zapicivaem srednee znachenie v T
    end
    T(i, 2)=mean(timeVector);
    K(i, 2)=k*ok;

%%%%%%%%%%%%%%%%%%%%% succesive over-relaxation %%%%%%%%%%%%%%%%%%%%%

    A = D{i}{1};
    b = D{i}{2};
    timeVector=zeros(N, 1);
    
    for j = 1 : 1 : N
        tic
    [x, ok, k] = my_successive_over_relaxation(A, b, x0, kmax, epsilon);
    timeVector(j, 1)=toc*ok;
    if ~ok   
        break
    end   
    end
    T(i, 3)=mean(timeVector);
    K(i, 3)=k*ok;
end







for i = 1 : 1 : 4
A = D{i}{1};
b = D{i}{2};
fprintf('A    ');
fprintf('\r\n');
[m, n]=size(A);
fprintf('razmernost A -- %d x %d', m, n);
fprintf('\r\n');
disp(A);
fprintf('b    ');
fprintf('\r\n');
n=length(b);
fprintf('razmersnost b -- %d x 1', n);
fprintf('\r\n');
disp(b);
%methods('Jacobi', 'Gauss_Seidel', 'Successive over-relaxation');
    Jacobi=X(:,1,i);
    Gauss_Seidel=X(:,2,i);
    Succ_over_relax=X(:,3,i);

    Tablica=table(Jacobi, Gauss_Seidel, Succ_over_relax);
disp(Tablica);
end


methods={'Jacobi', 'Gauss-Seidel', 'Successive over-relaxation'};
K=K';

subplot(1,2,1);
bar(T');
ax = gca;
ax.XTickLabels=methods;
title('T, time');
names={'A>0, Symmetric','A<0, Symmetric','A<0, Sparse','A<0, Full Sparse'};
length(D);
legend(names);

subplot(1,2,2);
bar(K);
ax = gca;
ax.XTickLabels=methods;
title('Iterations quantity');
legend(names);





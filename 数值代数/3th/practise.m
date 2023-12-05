%% 计算Householder向量 详情见House.m
clear all
%验证
x=[1;2;3;4];
n=size(x,1);
[beta,v]=House(x);
H=eye(n)-beta*(v*conj(v'));
H*x

%% Givens变换 详情见givens.m
clear all
%验证
x=[1;5];
[c,s]=givens(x(1),x(2));
G=[c,s;-s,c];
G*x

%% Gram-Schmidt过程
clear all
A=[1,2;5,6;9,10];%此处对A的要求是列满秩的
[m,n]=size(A);
q=zeros(m,n);r=zeros(n);
r(1,1)=norm(A(:,1));%norm函数在未添加其他参数是求的是二范数
q(:,1)=A(:,1)/r(1,1);
for j=2:n
    q(:,j)=A(:,j);
    for i=1:j-1
        r(i,j)=q(:,i)'*A(:,j);
        q(:,j)=q(:,j)-r(i,j)*q(:,i);
    end
    r(j,j)=norm(q(:,j));
    q(:,j)=q(:,j)/r(j,j);
end
q  %应为一个列正交矩阵
r  %应为一上三角矩阵

%% 基于MGS的QR分解
clear all
A=[1,2,2,4;5,6,6,5;9,10,10,5];
[m,n]=size(A);
R=zeros(n);
Q=zeros(m,n);
if A(:,1)~=0
    R(1,1)=norm(A(:,1));
    Q(:,1)=A(:,1)/norm(A(:,1));
end
for k=2:n
    Q(:,k)=A(:,k);
    for i=1:k-1
        R(i,k)=Q(:,i)'*Q(:,k);
        Q(:,k)=Q(:,k)-R(i,k)*Q(:,i);
    end
    if Q(:,k)~=0
        R(k,k)=norm(Q(:,k));
        Q(:,k)=Q(:,k)/R(k,k);
    end
end
Q
R

%% 基于Householder变换的QR分解
clear all
A=[1,2,4;5,6,5;9,10,5];%这里最好使m>=n,否则会使下面的for循环中的的x在某一次循环时为空值，从而报错
[m,n]=size(A);
% 可以在此处加入n=min(p[n,m])；以解决上面问题
Q=eye(m);
for k=1:n
    x=A(k:m,k);
    [beta,v]=House(x);
    A(k:m,k:n)=(eye(m-k+1)-beta*(v*v'))*A(k:m,k:n);
    Q(:,k:m)=Q(:,k:m)*(eye(m-k+1)-beta*(v*v'));
end
Q
A

%% 基于Givens变换的QR分解
clear all
A=[1,2,4;5,6,5;9,10,5];%m>=n
[m,n]=size(A);
Q=eye(m);
for k=1:n
    for i=k+1:m
        [c,s]=givens(A(k,k),A(i,k));
        G=[c,s;-s,c];
        GA=G*[A(k,k:n);A(i,k:n)];
        A(k,k:n)=GA(1,:);
        A(i,k:n)=GA(2,:);
        QGt=[Q(1:m,k),Q(1:m,i)]*G';
        Q(1:m,k)=QGt(:,1);
        Q(1:m,i)=QGt(:,2);
    end
end
Q
A

%% 分别用三种方法求解最小二乘问题，比较运算时间
clear all
close all
NN=500:500:5000;
length_NN=length(NN);
tt=zeros(length_NN,3);
for k=1:length_NN
    n=NN(k);
    m=2*n;
    A=randn(m,n);
    b=randn(m,1);
    fprintf('*******n=%d*******\n',n);
    %normal equation
    t0=clock;
    x1=(A'*A)\(A'*b);
    tt(k,1)=etime(clock,t0);
    
    % QR
    t0=clock;
    [Q,R]=qr(A,0);%[Q,R] = qr(A,0) 生成精简分解
    x2=R\(Q'*b);
    tt(k,2)=etime(clock,t0);
    fprintf('QR:relerr=%.2e\n',norm(x2-x1)/norm(x1));
    
    % SVD
    t0=clock;
    [U,S,V]=svd(A,'econ');
    x3=V*(diag(S).\(U'*b));
    tt(k,3)=etime(clock,t0);
    fprintf('SVD:relerr=%.2e\n\n',norm(x3-x1)/norm(x1));
end






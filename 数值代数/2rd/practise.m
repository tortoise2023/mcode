%% LU分解（KIJ型）
clear all
A=[1 2 3 4;5 6 7 8;9 10 12 12;13 14 15 8];
n=size(A,1);
L=eye(n);U=zeros(n);
for k=1:n-1
    for i=k+1:n
        L(i,k)=A(i,k)/A(k,k);%计算L的第k列
    end
    for j=k:n
        U(k,j)=A(k,j);%计算U的第k行
    end
    for i=k+1:n
        for j=k+1:n
            A(i,j)=A(i,j)-L(i,k)*U(k,j);
        end
    end
end
U(n,n)=A(n,n);
L
U

%% LU分解（用A存储L和U）（KIJ型）(mylu.m)
clear all
A=[1,2,3,4;5,6,7,8;9,10,16,12;13,14,15,17];
n=size(A,1);
for k=1:n-1
    for i=k+1:n
        A(i,k)=A(i,k)/A(k,k);
        for j=k+1:n
            A(i,j)=A(i,j)-A(i,k)*A(k,j);
        end
    end
end
A

%% LU分解（IKJ型）
clear all
A=[1 2 3 4;5 6 7 8;9 10 12 12;13 14 15 8];
n=size(A,1);
for i=2:n
    for k=1:i-1
        A(i,k)=A(i,k)/A(k,k);
        for j=k+1:n
            A(i,j)=A(i,j)-A(i,k)*A(k,j);
        end
    end
end
A

%% LU分解（待定系数法或Doolittle法）详情见mylu2.m

%% 向前回代求解Ly=b（行存储方式）
b=ones(n,1);
y=zeros(n,1);
y(1)=b(1)/L(1,1);
for i=2:n
    for j=1:i-1
        b(i)=b(i)-L(i,j)*y(j);
    end
    y(i)=b(i)/L(i,i);
end
y

%% 向前回代求解Ux=y（列存储方式）
x=zeros(n,1);
x(n)=y(n)/U(n,n);
for k=n:-1:1
    x(k)=y(k)/U(k,k);
    for i=k-1:-1:1
        y(i)=y(i)-x(k)*U(i,k);
    end
end
x

%% 部分选主元LU分解 详情见PLU.m
clear all
A=[1 2 3 4;5 6 7 8;9 10 12 12;13 14 15 8];
n=size(A,1);
L=zeros(n);U=zeros(n);
[A,p]=PLU(A);%
for i=1:n-1
    L=L+diag(diag(A,-i),-i);
    U=U+diag(diag(A,i),i);
end
L=L+eye(n)
U=U+diag(diag(A))
    
%% Cholesky分解（要求矩阵A对称正定）
clear all
A=[7 2 1 4;2 6 6 3;1 6 12 2;4 3 2 14];
n=size(A,1);
L=zeros(n);
for j=1:n
    L(j,j)=sqrt(A(j,j)-sum(L(j,1:j-1).^2));
    for i=j+1:n
        L(i,j)=(A(i,j)-sum(L(i,1:j-1)*L(j,1:j-1)'))/L(j,j);
    end
end
L

%% 改进的平方根法
%先计算分解
clear all
A=[7 2 1 4;2 6 6 3;1 6 12 2;4 3 2 14];
n=size(A,1);
L=zeros(n);d=zeros(n,1);
d(1)=A(1,1);
for i=2:n
    L(i,1)=A(i,1)/d(1);
end
for j=2:n
    d(j)=A(j,j)-L(j,1:j-1).^2*d(1:j-1);
    for i=j+1:n
        L(i,j)=(A(i,j)-L(i,1:j-1).*d(1:j-1)'*L(j,1:j-1)')/d(j);
    end
end
L=L+eye(n)
d
%解方程组：Ly=b 和 DL'x=y
b=ones(n,1);
y=zeros(n,1);
x=zeros(n,1);
y(1)=b(1);
for i=2:n
    y(i)=b(i)-L(i,1:i-1)*y(1:i-1);
end
x(n)=y(n)/d(n);
for i=n-1:-1:1
    x(i)=y(i)/d(i)-L(i+1:n,i)'*x(i+1:n);
end
x

%% 追赶法
clear all
A=[1 2 3 4;5 6 7 8;9 10 12 12;13 14 15 8];
n=size(A,1);
a=diag(A,-1);
b=diag(A);
c=diag(A,1);%三对角矩阵取自A
x=zeros(n,1);y=zeros(n,1);f=ones(n,1);%f为全为1的列向量
bata=zeros(n,1);alpha=zeros(n,1);
alpha(1)=a(1);
bata(1)=c(1)/b(1);
y(1)=f(1);
for i=2:n-1
    alpha(i)=b(i)-a(i-1)*bata(i-1);
    bata(i)=c(i)/alpha(i);
    y(i)=(f(i)-a(i-1)*y(i-1))/alpha(i);
end
alpha(n)=b(n)-a(n-1)*bata(n-1);
y(n)=(f(n)-a(n-1)*y(n-1))/alpha(n);
x(n)=y(n);
for i=n-1:-1:1
    x(i)=y(i)-bata(i)*x(i+1);
end
x
        
%% 求解Yule-Walker方程组的Durbin   
clear all
t=[2;3;4;5];
n=size(t,1);
x=zeros(n,1);x(1)=-t(1);bata=1;alpha=-t(1);
for k=1:n-1
    bata=(1-alpha^2)*bata;
    alpha=-(t(k+1)+t(1:k)'*x(k:-1:1))/bata;
    x(1:k)=x(1:k)+alpha.*x(k:-1:1);
    x(k+1)=alpha;
end
x
Tn=zeros(n);
for k=1:n-1
    Tn=Tn+diag(t(k)*ones(n-k,1),k);
end
Tn=Tn+eye(n)+Tn'
left=Tn*x
right=-t
    
%% 求解对称正定Toeplitz线性方程组的Levinson算法
clear all
t=[12;4;5];n=size(t,1)+1;
b=ones(n,1);y=zeros(n,1);x=zeros(n,1);
y(1)=-t(1);x(1)=b(1);bata=1;alpha=-t(1);
for k=1:n-1
    bata=(1-alpha^2)*bata;
    mu=(b(k+1)-t(1:k)'*x(k:-1:1))/bata;
    x(1:k)=x(1:k)+mu*y(k:-1:1);
    x(k+1)=mu;
    if(k<n-1)
        alpha=-(t(k+1)+t(1:k)'*y(k:-1:1))/bata;
        y(1:k)=y(1:k)+alpha*y(k:-1:1);
        y(k+1)=alpha;
    end
end
x
Tn=zeros(n);%构造对称正定Toeplitz矩阵进行验证
for k=1:n-1
    Tn=Tn+diag(t(k)*ones(n-k,1),k);
end
Tn=Tn+eye(n)+Tn';
left=Tn*x
right=b









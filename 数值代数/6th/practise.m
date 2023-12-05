%% Jacobi迭代
clear all
A=rand(10)+9*eye(10);
n=size(A,1);
b=ones(n,1);tol=1e-4;
x0=zeros(n,1);x0(1)=1;delta=norm(x0,2);
x1=zeros(n,1);k=0;
while delta>tol
    for i=1:n
        x1(i)=(b(i)-A(i,:)*x0+x0(i)*A(i,i))/A(i,i);
    end
    delta=norm(x0-x1,2);
    x0=x1;k=k+1;
end
x1

%% Gauss-Seidel迭代
clear all
A=rand(10)+10*eye(10);
n=size(A,1);
b=ones(n,1);tol=1e-4;
x0=zeros(n,1);x0(1)=1;delta=norm(x0,2);
x1=zeros(n,1);k=0;
while delta>tol
    for i=1:n
        x1(i)=1/A(i,i)*(b(i)-A(i,1:(i-1))*x1(1:(i-1))-A(i,(i+1):n)*x0((i+1):n));
    end
    delta=norm(x0-x1,2);
    x0=x1;k=k+1;
end
x1

%% 求解线性方程组的SOR迭代方法
clear all
A=rand(100)+100*eye(100);
n=size(A,1);
b=ones(n,1);tol=1e-4;
x0=zeros(n,1);x0(1)=1;delta=norm(x0,2);
x1=zeros(n,1);k=0;
omiga=1.5;
while delta>tol
    for i=1:n
        x1(i)=(1-omiga)*x0(i)+omiga/A(i,i)*(b(i)-A(i,1:(i-1))*x1(1:(i-1))-A(i,(i+1):n)*x0((i+1):n));
    end
    delta=norm(x0-x1,2);
    x0=x1;k=k+1;
end
x1

%% SSOR方法
clear all
A=rand(10)+10*eye(10);
n=size(A,1);
b=ones(n,1);tol=1e-4;
x0=zeros(n,1);x0(1)=1;delta=norm(x0,2);
x5=zeros(n,1);
x1=zeros(n,1);k=0;
omiga=1.5;
while delta>tol
    for i=1:n
        x5(i)=(1-omiga)*x0(i)+omiga/A(i,i)*(b(i)-A(i,1:(i-1))*x5(1:(i-1))-A(i,(i+1):n)*x0((i+1):n));
    end
    for i=n:-1:1
        x1(i)=(1-omiga)*x5(i)+omiga/A(i,i)*(b(i)-A(i,1:(i-1))*x5(1:(i-1))-A(i,(i+1):n)*x1((i+1):n));
    end
    delta=norm(x0-x1,2);
    x0=x1;k=k+1;
end
x1

%% 求解二维离散Poisson方程的Jacobi迭代方法
clear all
h=1e-2;
xlim=[0,1];ylim=[0,1];
nxstep=floor((xlim(2)-xlim(1))/h);nystep=floor((ylim(2)-ylim(1))/h);
u0=zeros(nystep+1,nxstep+1);
for i=1:nystep
    u0(i+1,nxstep+1)=G(xlim(2),(i-1)*h);
end
for i=1:nxstep
    u0(nystep+1,i)=G((i-1)*h,ylim(2));
    
end
u1=u0;
tol=1e-4;delta=1;k=0;
while delta>tol
    for i=2:nystep
        for j=2:nxstep
            u1(i,j)=1/4*(h^2*F(h*(j-1),h*(i-1))+u0(i+1,j)+u0(i-1,j)+u0(i,j+1)+u0(i,j-1));
        end
    end
    delta=norm(u0-u1,1);
    u0=u1;
    k=k+1;
end
x=xlim(1):h:xlim(2);y=ylim(1):h:ylim(2);
[X,Y]=meshgrid(x,y);
mesh(X,Y,u1)

%% 求解二维离散Possion方程的红黑排序G-S迭代方法
clear all
h=1e-2;
xlim=[0,1];ylim=[0,1];
nxstep=floor((xlim(2)-xlim(1))/h);nystep=floor((ylim(2)-ylim(1))/h);
u0=zeros(nystep+1,nxstep+1);
for i=1:nystep
    u0(i+1,nxstep+1)=G(xlim(2),(i-1)*h);
end
for i=1:nxstep
    u0(nystep+1,i)=G((i-1)*h,ylim(2));
end
u1=u0;
tol=1e-4;delta=1;k=0;
while delta>tol
    for i=2:nystep
        for j=2+mod(i,2):2:nxstep
            u1(i,j)=1/4*(h^2*F(h*(j-1),h*(i-1))+u0(i+1,j)+u0(i-1,j)+u0(i,j+1)+u0(i,j-1));
        end
    end
    for i=2:nystep
        for j=3+mod(i,2):2:nxstep
            u1(i,j)=1/4*(h^2*F(h*(j-1),h*(i-1))+u0(i+1,j)+u0(i-1,j)+u0(i,j+1)+u0(i,j-1));
        end
    end
    delta=norm(u0-u1,1);
    u0=u1;
    k=k+1;
end
x=xlim(1):h:xlim(2);y=ylim(1):h:ylim(2);
[X,Y]=meshgrid(x,y);
mesh(X,Y,u1)

%% 求解二维离散Possion的SOR迭代方法
clear all
h=1e-2;
xlim=[0,1];ylim=[0,1];
nxstep=floor((xlim(2)-xlim(1))/h);nystep=floor((ylim(2)-ylim(1))/h);
u0=zeros(nystep+1,nxstep+1);
for i=1:nystep
    u0(i+1,nxstep+1)=G(xlim(2),(i-1)*h);
end
for i=1:nxstep
    u0(nystep+1,i)=G((i-1)*h,ylim(2));
end
u1=u0;
tol=1e-4;delta=1;k=0;omiga=1.5;
while delta>tol
    for i=2:nystep
        for j=2:nxstep
            u1(i,j)=(1-omiga)*u0(i,j)+omiga/4*(h^2*F(h*(j-1),h*(i-1))+u0(i+1,j)+u1(i-1,j)+u0(i,j+1)+u1(i,j-1));
        end
    end
    delta=norm(u0-u1,1);
    u0=u1;
    k=k+1;
end
x=xlim(1):h:xlim(2);y=ylim(1):h:ylim(2);
[X,Y]=meshgrid(x,y);
mesh(X,Y,u0)

%% 求解二维离散Possion方程的红黑排序SOR迭代方法
clear all
h=1e-2;
xlim=[0,1];ylim=[0,1];
nxstep=floor((xlim(2)-xlim(1))/h);nystep=floor((ylim(2)-ylim(1))/h);
u0=zeros(nystep+1,nxstep+1);
for i=1:nystep
    u0(i+1,nxstep+1)=G(xlim(2),(i-1)*h);
end
for i=1:nxstep
    u0(nystep+1,i)=G((i-1)*h,ylim(2));
end
u1=u0;
tol=1e-4;delta=1;k=0;omiga=1.3;
while delta>tol
    for i=2:nystep
        for j=2+mod(i,2):2:nxstep
            u1(i,j)=(1-omiga)*u0(i,j)+omiga/4*(h^2*F(h*(j-1),h*(i-1))+u0(i+1,j)+u1(i-1,j)+u0(i,j+1)+u1(i,j-1));
        end
    end
    for i=2:nystep
        for j=3+mod(i,2):2:nxstep
            u1(i,j)=(1-omiga)*u0(i,j)+omiga/4*(h^2*F(h*(j-1),h*(i-1))+u0(i+1,j)+u1(i-1,j)+u0(i,j+1)+u1(i,j-1));
        end
    end
    delta=norm(u0-u1,1);
    u0=u1;
    k=k+1;
end
x=xlim(1):h:xlim(2);y=ylim(1):h:ylim(2);
[X,Y]=meshgrid(x,y);
mesh(X,Y,u1)
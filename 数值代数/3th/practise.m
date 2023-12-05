%% ����Householder���� �����House.m
clear all
%��֤
x=[1;2;3;4];
n=size(x,1);
[beta,v]=House(x);
H=eye(n)-beta*(v*conj(v'));
H*x

%% Givens�任 �����givens.m
clear all
%��֤
x=[1;5];
[c,s]=givens(x(1),x(2));
G=[c,s;-s,c];
G*x

%% Gram-Schmidt����
clear all
A=[1,2;5,6;9,10];%�˴���A��Ҫ���������ȵ�
[m,n]=size(A);
q=zeros(m,n);r=zeros(n);
r(1,1)=norm(A(:,1));%norm������δ�����������������Ƕ�����
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
q  %ӦΪһ������������
r  %ӦΪһ�����Ǿ���

%% ����MGS��QR�ֽ�
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

%% ����Householder�任��QR�ֽ�
clear all
A=[1,2,4;5,6,5;9,10,5];%�������ʹm>=n,�����ʹ�����forѭ���еĵ�x��ĳһ��ѭ��ʱΪ��ֵ���Ӷ�����
[m,n]=size(A);
% �����ڴ˴�����n=min(p[n,m])���Խ����������
Q=eye(m);
for k=1:n
    x=A(k:m,k);
    [beta,v]=House(x);
    A(k:m,k:n)=(eye(m-k+1)-beta*(v*v'))*A(k:m,k:n);
    Q(:,k:m)=Q(:,k:m)*(eye(m-k+1)-beta*(v*v'));
end
Q
A

%% ����Givens�任��QR�ֽ�
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

%% �ֱ������ַ��������С�������⣬�Ƚ�����ʱ��
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
    [Q,R]=qr(A,0);%[Q,R] = qr(A,0) ���ɾ���ֽ�
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






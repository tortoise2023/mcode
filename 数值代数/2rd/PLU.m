function [A,p]=PLU(A)
[~,n]=size(A);
p=1:n;%用于记录置换矩阵
for i=1:n-1
    [a,k]=max(abs(A(i:n,i)));%选列主元
    if a==0
        error('Error:第 %d 步的列主元为0！\n',i);
    end
    k=k+i-1;
    if k~=i
        tmp=A(i,:);A(i,:)=A(k,:);A(k,:)=tmp;%交换A的第i行与第k行
        tmp=p(i);p(i)=p(k);p(k)=tmp;%更新置换矩阵
    end
    A(i+1:n,i)=A(i+1:n,i)/A(i,i);%计算L的第i列
    A(i+1:n,i+1:n)=A(i+1:n,i+1:n)-A(i+1:n,i)*A(i,i+1:n);%更新A(i+1:n,i+1:n)`
end
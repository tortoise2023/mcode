function [A,p]=PLU(A)
[~,n]=size(A);
p=1:n;%���ڼ�¼�û�����
for i=1:n-1
    [a,k]=max(abs(A(i:n,i)));%ѡ����Ԫ
    if a==0
        error('Error:�� %d ��������ԪΪ0��\n',i);
    end
    k=k+i-1;
    if k~=i
        tmp=A(i,:);A(i,:)=A(k,:);A(k,:)=tmp;%����A�ĵ�i�����k��
        tmp=p(i);p(i)=p(k);p(k)=tmp;%�����û�����
    end
    A(i+1:n,i)=A(i+1:n,i)/A(i,i);%����L�ĵ�i��
    A(i+1:n,i+1:n)=A(i+1:n,i+1:n)-A(i+1:n,i)*A(i,i+1:n);%����A(i+1:n,i+1:n)`
end
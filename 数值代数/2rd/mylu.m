function A=mylu(A)
n=size(A,1);%LU·Ö½â£¨KIJÐÍ£©
for k=1:n-1
    if A(k,k)==0
        fprintf('Error:A(%d,%d)=0!\n',k,k);
        return;
    end
    A(k+1:n,k)=A(k+1:n,k)/A(k,k);
    A(k+1:n,k+1:n)=A(k+1:n,k+1:n)-A(k+1:n,k)*A(k,k+1:n);
end
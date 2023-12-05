function A=mylu2(A)
[~,n]=size(A);
for k=1:n
    A(k,k)=A(k,k)-A(k,1:k-1)*A(1:k-1,k);
    if (A(k,k)==0)
        fprintf('Error:A(%d,%d)=0!\n',k,k);return;
    end
    A(k,k+1:n)=A(k,k+1:n)-A(k,1:k-1)*A(1:k-1,k+1:n);
    A(k+1:n,k)=(A(k+1:n,k)-A(k+1:n,1:k-1)*A(1:k-1,k))/A(k,k);
end
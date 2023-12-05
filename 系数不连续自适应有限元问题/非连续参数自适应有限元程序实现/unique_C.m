function L=unique_C(line)
n=size(line,2);
L=zeros(size(line,1),n);
L(:,1)=line(:,1);
L_num=1;
for i=2:n
    k=1;
    for j=L_num:-1:1
        if line(:,i)==L(:,j)
            k=0;
            break
        end
    end
    if k==1
        L_num=L_num+1;
        L(:,L_num)=line(:,i);
    end
end
L=L(:,1:L_num);




function V=Vk(p,t,a,K)
center_tri=t(1:3,K);
flag=0;
for i=1:3
    if p(1,t(i,K))==0 && p(2,t(i,K))==0
        flag=1;
    end
end

if flag==1
    

    all_tri=[];
    for i=1:3
        all_tri=[all_tri;find(t(1:3,:)==center_tri(i))];
    end
    all_tri=mod(all_tri-1,size(t,2))+1;
    all_tri=unique_C(all_tri');
    %计算所有的三角形上的值
    V_all=zeros(1,size(all_tri,2));
    for i=1:size(all_tri,2)
        
        
        mean_tri=sum(p(:,t(1:3,all_tri(i))'),2)/3;
        V_all(i)=a(mean_tri(1),mean_tri(2));
    end
    
    center_mean=sum(p(:,t(1:3,K)'),2)/3;
    a_center=a(center_mean(1),center_mean(2));
    V=a_center/min(V_all);
else
    V=1;
end



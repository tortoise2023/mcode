function longline=TheLong(p,t,K)
line=[t(1:2,K)',0;t(2:3,K)',0;t(3,K),t(1,K),0];
for i=1:3
    line(i,3)=norm(p(:,line(i,1))-p(:,line(i,2)));
end
line=sortrows(line,3,'descend');
longline=line(1,:);
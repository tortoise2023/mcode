function DrawExactSolution()
gamma=0.1;
[~,rou,sigma]=KelloggParameter(gamma);
u_bound=@(x,y)u_r_theta(x,y,rou,sigma,gamma);


h=0.01;
x=-1:h:1;
y=-1:h:1;

z=zeros(1,size(y,2));
for i=1:size(x,2)

    for k=1:size(y,2)
        z(k)=u_bound(x(i),y(k));
    end
    plot3(x(i)*ones(1,size(y,2)),y,z)
    hold on
    z=zeros(1,size(x,2));
end
for i=1:size(y,2)

    for k=1:size(x,2)
        z(k)=u_bound(x(k),y(i));
    end
    plot3(x,y(i)*ones(1,size(x,2)),z)
    z=zeros(1,size(x,2));

end


function AdaptiveDiscontinusCoefficientsSolver2D()
%clear all
gamma=0.1;
% R=5.82842712474619;
% rou=pi/4;
% sigma=-2.35619449019;
[R,rou,sigma]=KelloggParameter(gamma);
u_bound=@(x,y)u_r_theta(x,y,rou,sigma,gamma);

f=inline('0','x','y');
g=Rectg(-1,-1,1,1);
aa=@(x,y)a(x,y,R);
[p,e,t]=Generate_the_initial_grid();
theta=0.1;
ErrorAndFreedom=[];
%for l=1:2

while size(p,2)<5000
    b=LoadAssembler2D(p,t,f);
    A=StiffnessAssembler2D(p,t,aa);
    %加入边界条件
    out_nodes=e(1,:);
    for i=1:size(out_nodes,2)
        A(out_nodes(i),:)=0;
        A(out_nodes(i),out_nodes(i))=1;
        b(out_nodes(i))=u_bound(p(1,out_nodes(i)),p(2,out_nodes(i)));
    end
    

    xi=A\b;
    %计算与真解的误差
    xi_true=zeros(size(xi,1),1);
    for j=1:size(xi_true)
        xi_true(j)=u_bound(p(1,j),p(2,j));
    end
    %Error=xi_true-xi;
    
    
     
%     xlabel('x'),ylabel('y'),title('u-u_h',size(t,2))
%     figure(2),pdesurf(p,t,xi)
%     xlabel('x'),ylabel('y'),title('u_h',size(t,2))

    

    [elements,E]=LocalErrorIndicator(p,e,t,theta,aa,f,xi);

%     eta=pdejmps(p,t,1,0,1,xi,1,1,1);
%     tol=0.95*max(eta);
%     elements=find(eta>tol);
%     
    %EEE=max(Error);
    
    ErrorAndFreedom=[ErrorAndFreedom;E,size(xi,1)];
    %[p,e,t]=RefinemeshChange(p,e,t,elements);
    [p,e,t]=refinemesh(g,p,e,t,elements','regular');
    
%      figure(4),pdemesh(p,e,t)
%      title('number of point',size(p,2))
    %pause(0.5)
end


figure(3);hold on
plot(log(ErrorAndFreedom(:,2))/log(10),log(ErrorAndFreedom(:,1))/log(10))
ylabel('log_{10}(||u-u_h||_\infty)'),xlabel('log_{10}(DOFs)')

function [elements,E]=LocalErrorIndicator(p,e,t,theta,a,f,xi)
line=[t(1,:),t(2,:),t(3,:);t(2,:),t(3,:),t(1,:)];
line=[e(1:2,:),line];
line=sort(line);
line=unique_C(line);
inter_line=line(:,size(e,2)+1:end);
eta=zeros(1,size(inter_line,2));
connect_tri=[];
for i=1:size(inter_line,2)
    nn=p(:,inter_line(1,i))-p(:,inter_line(2,i));
    %He=norm(nn);
    
    

    %找第一个点符合的三角形
    w_tri=find(t(1:3,:)'==inter_line(1,i));
    w_tri=mod(w_tri-1,size(t,2))+1;
    %在第一个点符合的基础上找第二个点符合的三角形
    K=find(t(1:3,w_tri')'==inter_line(2,i));
    K=mod(K-1,size(w_tri,1))+1;
    K=w_tri(K);
    connect_tri=[connect_tri;K'];

    %计算第一个三角形上xi的梯度
    localtion_1=p(:,t(1:3,K(1))');
    A=[localtion_1',xi(t(1:3,K(1)))];
    A(2,:)=A(2,:)-A(1,:);
    A(3,:)=A(3,:)-A(1,:);
    x_1=[A(2,2)*A(3,3)-A(2,3)*A(3,2);A(2,3)*A(3,1)-A(2,1)*A(2,3)];
    x_1=-x_1/(A(2,1)*A(3,2)-A(2,2)*A(3,1));
    %计算eta的第一项
    B=zeros(2);
    B(1,:)=A(1,1:2)-A(3,1:2);
    B(2,:)=A(2,1:2)-A(2,1:2);
    S_1=abs(det(B))/2;
    Vk_1=Vk(p,t,a,K(1));
    tri_length_1=[norm(localtion_1(:,1)-localtion_1(:,2))
        norm(localtion_1(:,3)-localtion_1(:,2))
        norm(localtion_1(:,1)-localtion_1(:,3))];
    Hk_1=max(tri_length_1);
    localtion_mean_1=sum(localtion_1,2)./3;
    f_meam_1=f(localtion_mean_1(1),localtion_mean_1(2));
    a_1=a(localtion_mean_1(1),localtion_mean_1(2));
    eta_1=Vk_1*Hk_1/a_1*f_meam_1^2*S_1;

    

    %计算第二个三角形上xi的梯度
    localtion_2=p(:,t(1:3,K(2))');
    A=[localtion_2',xi(t(1:3,K(2)))];
    A(2,:)=A(2,:)-A(1,:);
    A(3,:)=A(3,:)-A(1,:);
    x_2=[A(2,2)*A(3,3)-A(2,3)*A(3,2);A(2,3)*A(3,1)-A(2,1)*A(2,3)];
    x_2=-x_2/(A(2,1)*A(3,2)-A(2,2)*A(3,1));
    %计算eta的第二项
    B=zeros(2);
    B(1,:)=A(1,1:2)-A(3,1:2);
    B(2,:)=A(2,1:2)-A(2,1:2);
    S_2=abs(det(B))/2;
    Vk_2=Vk(p,t,a,K(2));
    tri_length_2=[norm(localtion_2(:,1)-localtion_2(:,2))
        norm(localtion_2(:,3)-localtion_2(:,2))
        norm(localtion_2(:,1)-localtion_2(:,3))];
    Hk_2=max(tri_length_2);
    localtion_mean_2=sum(localtion_2,2)./3;
    f_meam_2=f(localtion_mean_2(1),localtion_mean_2(2));
    a_2=a(localtion_mean_2(1),localtion_mean_2(2));
    eta_2=Vk_2*Hk_2/a_2*f_meam_2^2*S_2;


    %计算法向量
    normal_vector=zeros(1,2);
    normal_vector(1)=nn(2);
    normal_vector(2)=-nn(1);
    %normal_vector=normal_vector;
    %计算跳量
    Je=abs(normal_vector*(x_1*a_1-x_2*a_2));
    %计算Ve
    Ve=max(Vk_1,Vk_2);
    %计算eta的第三项
    ae=max(a_1,a_2);
    eta_3=Je^2*Ve/ae;
    %合并
    eta(i)=eta_1+eta_2+eta_3;
end

all_elements=[eta',connect_tri];
all_elements=sortrows(all_elements,'descend');
refine_num=1;
just_ok_eta=sum(eta,2)*theta^2;
for i=2:size(eta,2)
    all_elements(i,1)=all_elements(i,1)+all_elements(i-1,1);
    if all_elements(i,1)>just_ok_eta
        
        break
    end
    refine_num=refine_num+1;
end
elements=all_elements(1:refine_num,2:3);
elements=reshape(elements,1,[]);
elements=unique_C(elements);
E=sqrt(sum(eta));



    









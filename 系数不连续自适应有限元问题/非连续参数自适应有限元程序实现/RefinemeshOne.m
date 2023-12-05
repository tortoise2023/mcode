function [p,e,t,finish_K]=RefinemeshOne(p,e,t,K,pp)
finish_K=K;

longline=TheLong(p,t,K);
%判断最长边是否为边界边
outline_1=find(e(1:2,:)'==longline(1));
outline_1=mod(outline_1-1,size(e,2))+1;
outline_2=find(e(1:2,outline_1)'==longline(2));
outline_2=outline_1(mod(outline_2-1,size(outline_1,1))+1);

%细化p，t
p=[p,(p(:,longline(1))+p(:,longline(2)))/2];
if pp==0
    t=[t,t(:,K)];
    t(t(1:3,K)==longline(1),K)=size(p,2);
    t(t(1:3,end)==longline(2),end)=size(p,2);
else
    t=[t,t(:,K),t(:,K)];
    if longline(2)==pp(2)
        t(1,K)=longline(1);
        t(2,K)=size(p,2);
        t(3,K)=pp(1);

        t(1,end-1)=size(p,2);
        t(2,end-1)=longline(2);
        t(3,end-1)=size(p,2)-1;

        t(1,end)=size(p,2);
        t(2,end)=size(p,2)-1;
        t(3,end)=pp(1);
    else
        t(1,K)=longline(2);
        t(2,K)=pp(2);
        t(3,K)=size(p,2);

        t(1,end-1)=size(p,2);
        t(2,end-1)=pp(2);
        t(3,end-1)=size(p,2)-1;

        t(1,end)=size(p,2);
        t(2,end)=size(p,2)-1;
        t(3,end)=pp(1);
       
    end
end



if size(outline_2,1)==1
    %细化e
    e=[e,e(:,outline_2)];
    e(2,outline_2)=size(p,2);
    e(4,outline_2)=(e(4,outline_2)+e(3,outline_2))/2;
    e(1,end)=size(p,2);
    e(3,end)=e(4,outline_2);
    
else
    %找到与最长边相邻的另一个三角形

    tri_num=find(t(1:3,:)'==longline(1));
    tri_num=mod(tri_num-1,size(t,2))+1;
    tri_num_true=find(t(1:3,tri_num)'==longline(2));
    tri_num_true=tri_num(mod(tri_num_true-1,size(tri_num,1))+1);

    if size(tri_num_true,1)==1

        %看是否另一个三角形的最长边也是这个
        
        longline_2=TheLong(p,t,tri_num_true);
        if longline_2(3)==longline(3)
            %细化另一个三角形的p，t
            t=[t,t(:,tri_num_true)];
            t(t(1:3,tri_num_true)==longline_2(1),tri_num_true)=size(p,2);
            t(t(1:3,end)==longline_2(2),end)=size(p,2);
            finish_K=[finish_K,tri_num_true];
        else
            [p,e,t,finish_K_K]=RefinemeshOne(p,e,t,tri_num_true,longline);
            finish_K=[finish_K,finish_K_K];
        end
    end
end
    



    



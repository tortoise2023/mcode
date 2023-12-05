function [p,e,t]=RefinemeshChange(p,e,t,elements)
while size(elements,2)>0
%for i=1:2
    
    [p,e,t,finish_K]=RefinemeshOne(p,e,t,elements(1),0);
    n=size(finish_K,2);
    elements=[finish_K,elements];
    elements=unique_C(elements);
    elements(1:n)=[];
end








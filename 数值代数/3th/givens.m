function [c,s]=givens(a,b)
if b==0
    if a>=0
        c=1;s=0;
    else
        c=-1;s=0;
    end
else
    if abs(b)>abs(a)
        tau=a/b;s=sign(b)/sqrt(1+tau^2);c=s*tau;
    else
        tau=b/a;c=sign(a)/sqrt(1+tau^2);s=c*tau;
    end
end
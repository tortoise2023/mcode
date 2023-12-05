function [R,rou,sigma]=KelloggParameter(gamma)
F=@(R,rou,sigma)[R+tan((pi/2-sigma)*gamma)*cot(rou*gamma);
    1/R+tan(rou*gamma)*cot(sigma*gamma);
    R+tan(sigma*gamma)*cot((pi/2-rou)*gamma)];
x0=[161.4476;pi/4;-14.9225];
X=NewtonN(F,x0,3);
R=X(1);
rou=X(2);
sigma=X(3);
%此参数使用牛顿迭代时极不可靠，猜测是由于函数不连续所造成的
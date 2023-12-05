function u=u_r_theta(x,y,rou,sigma,gamma)
r=sqrt(x^2+y^2);
if r==0
    u=0;
else
    if y>=0
        theta=acos(x/r);
    else
        theta=-acos(x/r)+2*pi;
    end
    if theta>=0 && theta<=pi/2
        mu=cos((pi/2-sigma)*gamma)*cos((theta-pi/2+rou)*gamma);
    elseif theta>=pi/2 && theta<=pi
        mu=cos(rou*gamma)*cos((theta-pi+sigma)*gamma);
    elseif theta>=pi && theta<=3*pi/2
        mu=cos(sigma*gamma)*cos((theta-pi-rou)*gamma);
    else
        mu=cos((pi/2-rou)*gamma)*cos((theta-3*pi/2-sigma)*gamma);
    end
    u=r^gamma*mu;
end


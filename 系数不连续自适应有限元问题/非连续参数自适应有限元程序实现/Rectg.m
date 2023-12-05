function r=Rectg(xmin,ymin,xmax,ymax)
r=[2,xmin,xmax,ymax,ymax,1,0;
    2,xmax,xmax,ymax,ymin,1,0;
    2,xmax,xmin,ymin,ymin,1,0;
    2,xmin,xmin,ymin,ymax,1,0]';

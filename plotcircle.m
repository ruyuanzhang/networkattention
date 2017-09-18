function plotcircle(x,y,r)
%x and y are the coordinates of the center of the plotcircle
%r is the radius of the plotcircle
%0.01 is the angle step, bigger values will draw the plotcircle faster but
%you might notice imperfections (not very smooth)
ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
plot(x+xp,y+yp);
end

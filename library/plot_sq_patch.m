function plot_sq_patch(x,y)

xcor = [x-hours(0.5) x+hours(0.5) x+hours(0.5) x-hours(0.5)];
ycor = [y(1) y(1) y(2) y(2)];
fill(xcor,ycor,rgb('coral'),'FaceAlpha',0.1,'LineStyle','none',...
     'Marker','none')

end
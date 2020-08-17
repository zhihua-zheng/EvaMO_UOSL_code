function plot_pdf(xdata1,xdata2,numB,xString)

colorString1 = 'azure';
colorString2 = 'coral';
off_ratio  = 0.03;

hisf1 = figure;
his1  = histogram(xdata1,numB,'Normalization','pdf');
pdf1  = his1.Values;
npdf1 = pdf1/max(pdf1);
xi1   = his1.BinEdges;
x1    = (xi1(1:end-1)+xi1(2:end))/2;

hisf2 = figure;
his2  = histogram(xdata2,numB,'Normalization','pdf');
pdf2  = his2.Values;
npdf2 = pdf2/max(pdf2);
xi2   = his2.BinEdges;
x2    = (xi2(1:end-1)+xi2(2:end))/2;

figure('position',[0 0 400 380]) 
plot(x1,npdf1,'color',rgb(colorString1),'linewidth',2.5)
hold on; grid on; ylim([0 1.15])
fill([x1 flip(x1)],[zeros(1,his1.NumBins) flip(npdf1)],rgb(colorString1),...
     'LineStyle','none','FaceAlpha',.25)
         
plot(x2,npdf2,'color',rgb(colorString2),'linewidth',2.5)
fill([x2 flip(x2)],[zeros(1,his2.NumBins) flip(npdf2)],rgb(colorString2),...
     'LineStyle','none','FaceAlpha',.25)

set(gca,'TickDir','out','fontsize',12)
ax = gca;
ax.Position = ax.Position + off_ratio*[1 1 -1 -1];
ylh = ylabel('Normalized PDF','fontsize',16);
xlh = xlabel(xString,'fontsize',16,'Interpreter','tex');
xlh.Units = 'Normalized';
ylh.Units = 'Normalized';
set(xlh,'Position',xlh.Position + [0 -off_ratio/6 0]);
set(ylh,'Position',ylh.Position + [-off_ratio/3 0 0]);

close(hisf1)
close(hisf2)

end
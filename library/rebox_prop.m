function rebox_prop(boxG,x,g,gv,uperc,cstr)

%% Find handles to boxplot components

uw  = findobj(boxG,'tag','Upper Whisker');        % handle to "Upper Whisker" line
uav = findobj(boxG,'tag','Upper Adjacent Value'); % handle to "Upper Adjacent Value" line
lw  = findobj(boxG,'tag','Lower Whisker');        % handle to "Lower Whisker" line
lav = findobj(boxG,'tag','Lower Adjacent Value'); % handle to "Lower Adjacent Value" line

otlh = findobj(boxG,'tag','Outliers'); % handle to "Outliers"
boxh = findobj(boxG,'tag','Box');      % handle to "Box"
medh = findobj(boxG,'tag','Median');   % handle to "Median"

%% Separate data into different groups

ng = numel(gv);
x_grouped  = cell(ng,1);
x_Uwhisker = cell(ng,1);
x_Lwhisker = cell(ng,1);
x_outlier  = cell(ng,1);
x_mean     = cell(ng,1);

lperc = 1 - uperc;

for k = 1:ng
    x_in_grp      = x(g == gv(ng-k+1));
    x_grouped{k}  = x_in_grp;
    x_mean{k}     = nanmean(x_grouped{k});
    x_Uwhisker{k} = quantile(x_grouped{k},uperc);
    x_Lwhisker{k} = quantile(x_grouped{k},lperc);
    x_outlier{k}  = x_in_grp(x_in_grp < x_Lwhisker{k} | ...
                             x_in_grp > x_Uwhisker{k});
end

%% Modify whisker lines, boxes, medians

for k = 1:ng
    
    uw(k).YData(1,2) = x_Uwhisker{k};
    uw(k).LineStyle  = '-';
    uw(k).LineWidth  = 1;
    
    uav(k).YData(:)  = x_Uwhisker{k};
    uav(k).LineWidth = 1.5;
    
    lw(k).YData(1,1) = x_Lwhisker{k};
    lw(k).LineStyle  = '-';
    lw(k).LineWidth  = 1;
    
    lav(k).YData(:)  = x_Lwhisker{k};
    lav(k).LineWidth = 1.5;
    
    boxh(k).LineWidth = 1;
    boxh(k).Color     = [0 0 0];
    patch(get(boxh(k),'XData'),get(boxh(k),'YData'),rgb(cstr),'FaceAlpha',.5);
    
    medh(k).LineWidth = 2;
    medh(k).Color     = [0 0 0];
end

%% Adjust outliers accordingly and add mean value

for k = 1:ng
    
    XData = unique(uw(k).XData);       
    plot(XData,x_mean{k},'d','MarkerFaceColor','none',...
        'MarkerEdgeColor','k','MarkerSize',5)
%     plot([min(get(boxh(k),'XData')) max(get(boxh(k),'XData'))],...
%          x_mean{k}*[1 1],'--','LineWidth',1,'color','k')

    delete(otlh(k))
    nout = numel(x_outlier{k});
    if nout == 0
        otlh(k) = plot(XData,NaN);
    else
        otlh(k) = plot(XData*ones(1,nout),x_outlier{k},'+',...
                  'color',rgb(cstr),'MarkerSize',4);
    end
end

end
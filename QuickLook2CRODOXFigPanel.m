function [outpt,k,kerr]=QuickLook2CRODOXFigPanel
% this one is used to make the figure


fils=dir('*WT_18*mat');

ftmp=load(fils(1).name);
out(1).D=ftmp.D;
%data=reshape(ftmp.data,4,12,[]);
%data=squeeze(mean(data));  %load data and average every 4 rows together
out(1).data=ftmp.data;
out(1).file=fils(1).name;

fils=dir('*day*18hrs*mat');

for i=1:length(fils)
   
    ftmp=load(fils(i).name);
    out(i+1).D=ftmp.D;
%     data=reshape(ftmp.data,4,12,[]);
% data=squeeze(mean(data));  %load data and average every 4 rows together
    out(i+1).data=ftmp.data;
    out(i+1).file=fils(i).name;
    
end


wtmean=mean(out(1).data);
wtmean=wtmean/wtmean(1);


ct=1;
for i=1:length(out)
   
    for j=1:length(out(i).data)
       
        [ccout,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(out(i).D,out(i).data(j,:)/out(i).data(j,1),@(cc,xx) 1./(1+(xx./cc(1)).^cc(2)),[100,2]);
        ci = nlparci(ccout,R,'covar',CovB);
        k(i,j)=ccout(1);
        klow(i,j)=ci(1,1);
        khigh(i,j)=ci(1,2);
% %       subplot(4,12,ct)
% %       if i>1
% %       plot(out(i-1).D,out(i-1).data(j,:)/out(i-1).data(j,1),'-x')
% %       hold on
% %       plot(out(1).D,wtmean,'r')
% %      ct=ct+1;   
% %      set(gca,'XTick',[])
% %      set(gca,'YTick',[])
% %      xlim([0 650])
% %      ylim([0.1 1])
% %       end
    end
    
end

k(abs(imag(k))>1e-10)=NaN;
% let's average over the ic50's for the 3 technical reps
for i=1:size(k,1)
    ktmp=k(i,:);  %48 ic50 vals
    ksq=reshape(ktmp,4,[]);
    kmean=mean(ksq);
    kMn(i,:)=kmean;
    kStd(i,:)=std(ksq)./sqrt(4);
      ksqcube{i}=ksq;
    
end

k=kMn;
kerr=kStd;

% figure
% subplot(1,2,1)
% for i=1:length(out)
%    plot(out(i).D,mean(out(i).data)/mean(out(i).data(:,1)),'-','linewidth',2)
%     hold on
%     
% % end
% 
% axis square
% set(gca,'fontsize',12)
% xlabel('[DOX] (\mug/mL)')
% ylabel('Cell Density')
% legend('WT','D2','D4','D6','D8')
% title('CRO-Selected Isolates')
% set(gca,'XScale','log')
% xlim([0 650])

%subplot(1,2,2)
% % % % % for i=1:5
% % % % %   h=errorbar(i+randn(size((k(i,:))))*0.03,k(i,:),kerr(i,:),'*','markersize',2,'linewidth',0.5);
% % % % %   alpha = 0.05;   
% % % % % % Set transparency (undocumented)
% % % % % set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha])
% % % % % 
% % % % %    hold on
% % % % %    errorbar(i,mean(k(i,:)),std(k(i,:))/sqrt(length(k(i,:))),'sk','markersize',10,'markerfacecolor','k')
% % % % %    plot(i,max(k(i,:)),'sk','markersize',10,'markerfacecolor','w')
% % % % % end
% % % % % axis square
% % % % % set(gca,'xtick',[1:6])
% % % % % set(gca,'XTickLabel',{'WT','D2','D4','D6','D8'})
% % % % %     set(gca,'fontsize',12)
% % % % %     ylabel('Estimated IC_{50}')
% % % % % plot([0,6],[1 1]*mean(k(1,:)),'--k','linewidth',1)
% % % % % 
% % % % % 
% % % % % for i=1:4
% % % % %    [hval,pval]=ttest2(k(1,:),k(i+1,:),'Vartype','unequal')
% % % % % end
    
wtMean=nanmean(k(1,:));
wtStErr=nanstd(k(1,:))/sqrt(sum(k(1,:)>0));
shadedErrorBar([0 6],[1 1],1*[1 1]*wtStErr/wtMean,'patchsaturation',0.01,'lineProps',{'k-','linestyle','none'})
hold on
for i=1:5
    
    hold on
  h=errorbar(i+randn(size((k(i,:))))*0.1,k(i,:)/wtMean,kerr(i,:)/wtMean,'o','markersize',4);
%   alpha = 0.05;   
% % Set transparency (undocumented)
% set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha])

   hold on
   errorbar(i,nanmean(k(i,:))/wtMean,nanstd(k(i,:))/sqrt(length(k(i,:)))/wtMean,'sk','markersize',10,'markerfacecolor','w','linewidth',2)
   %plot(i,nanmax(k(i,:)),'sk','markersize',10,'markerfacecolor','w')
   
end
axis square
set(gca,'xtick',[1:6])
set(gca,'XTickLabel',{'WT','2','4','6','8'})

    set(gca,'fontsize',12)
    ylabel('Estimated IC_{50}')
%plot([0,6],[1 1]*nanmean(k(1,:)),'--k','linewidth',1)

xlim([1.5 5.5])
xlabel('Days')
set(gca,'fontsize',12)


disp('CRO-DOX')
% [hval1stLastTime,pval1stLastTime]=ttest2(k(2,:),k(5,:),'Vartype','unequal')

  [pp,tbl,stats]=anova1(k(2:5,:)');

[results, means, h, nms]=multcompare(stats);
tbl = array2table(results,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
clusters=statGroupClusters(results,means) 
   ppAnova=pp

   writetable(tbl,'AnovaMultCompCRODOX.csv')
% for i=1:4
%    [hval,pval]=ttest2(k(1,:),k(i+1,:),'Vartype','unequal')
% end
    

set(gca,'YScale','log')
ylim([6e-1 1.6])
title('CRO-DOX')

outpt.wtMean=wtMean;
outpt.wtStErr=wtStErr;
outpt.k=k;
outpt.kerr=kerr;
outpt.kfull=ksqcube;
end
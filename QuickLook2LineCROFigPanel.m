function [outpt,k,kerr]=QuickLook2LineCROFigPanel
% this version is to make the figure panel

fils=dir('*WT_16*mat');

ftmp=load(fils(1).name);
out(1).D=ftmp.D;
%data=reshape(ftmp.data,4,12,[]);
%data=squeeze(mean(data));  %load data and average every 4 rows together
out(1).data=ftmp.data;
out(1).file=fils(1).name;

fils=dir('*day*16hrs*mat');

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

%figure
ct=1;

for i=1:length(out)
   
    for j=1:length(out(i).data)
        xtmp=out(i).D;
        ytmp=out(i).data(j,:)/out(i).data(j,1) ;   
        try
        [ccout]=nlinfit(xtmp,ytmp,@(cc,xx) 1./(1+(xx./cc(1)).^cc(2)),[300,2]);
       % ci = nlparci(ccout,R,'covar',CovB);
        k(i,j)=ccout(1);
        catch
            k(i,j)=NaN;
            disp('fit failed')
        end
%         klow(i,j)=ci(1,1);
%         khigh(i,j)=ci(1,2);
%       subplot(5,12,ct)
%    
%       plot(out(i).D,out(i).data(j,:)/out(i).data(j,1),'*-','linewidth',2)
%       hold on
%       plot(out(1).D,wtmean,'k:','LineWidth',2)
%      
%       if mod(j,4)==0
%       ct=ct+1;   
%       end
%      set(gca,'XTick',[])
%      set(gca,'YTick',[])
%      xlim([0 650])
%      ylim([0.1 1])
%    
    end
    
end

k(abs(imag(k))>1e-5)=NaN;

% let's average over the ic50's for the 3 technical reps
for i=1:size(k,1)
    ktmp=k(i,:);  %48 ic50 vals
    ksq=reshape(ktmp,4,[]);
    kmean=nanmean(ksq);
    kMn(i,:)=kmean;
    kStd(i,:)=nanstd(ksq)./sqrt(4);
    ksqcube{i}=ksq;
      
end

k=kMn;
kerr=kStd;

% % % figure
% % % subplot(1,2,1)
% % % for i=1:length(out)
% % %    plot(out(i).D,nanmean(out(i).data)/nanmean(out(i).data(:,1)),'-','linewidth',2)
% % %     hold on
% % %     
% % % end
% % % 
% % % axis square
% % % set(gca,'fontsize',12)
% % % xlabel('[CRO] (\mug/mL)')
% % % ylabel('Cell Density')
% % % legend('WT','D2','D4','D6','D8')
% % % title('LZD-Selected Isolates')
% % % set(gca,'XScale','log')
% % % xlim([0 650])
% % % 
% % % subplot(1,2,2)
% % % % for i=1:5
% % % %   h=errorbar(i+randn(size((k(i,:))))*0.03,k(i,:),kerr(i,:),'*','markersize',2,'linewidth',0.5);
% % % %   alpha = 0.15;   
% % % % % Set transparency (undocumented)
% % % % set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha])
% % % % 
% % % %    hold on
% % % %    errorbar(i,nanmean(k(i,:)),nanstd(k(i,:))/sqrt(length(k(i,:))),'sk','markersize',10,'markerfacecolor','k')
% % % % end
% % % % axis square
% % % % set(gca,'xtick',[1:6])
% % % % set(gca,'XTickLabel',{'WT','D2','D4','D6','D8'})
% % % %     set(gca,'fontsize',12)
% % % %     ylabel('Estimated IC_{50}')
% % % % plot([0,6],[1 1]*nanmean(k(1,:)),'--k','linewidth',1)
% % % % xlim([0.5 5.5])
% % % % 
% % % % for i=1:4
% % % %    [hval,pval]=ttest2(k(1,:),k(i+1,:),'Vartype','unequal')
% % % % end
% % % %     

wtMean=nanmean(k(1,:));
wtStErr=nanstd(k(1,:))/sqrt(sum(k(1,:)>0));
shadedErrorBar([0 6],[1 1],1*[1 1]*wtStErr/wtMean,'patchsaturation',0.01,'lineProps',{'k-','linestyle','none'})
hold on

outpt.wtMean=wtMean;
outpt.wtStErr=wtStErr;
outpt.k=k;
outpt.kerr=kerr;
outpt.kfull=ksqcube;

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

% for i=1:4
%    [hval,pval]=ttest2(k(1,:),k(i+1,:),'Vartype','unequal')
% end
%     

disp('LZD-CRO')
% [hval1stLastTime,pval1stLastTime]=ttest2(k(2,:),k(5,:),'Vartype','unequal')
 
 kt=k;
 kt(imag(k)>1e-2)=NaN;
 kt=real(kt);
 [pp,tbl,stats]=anova1((kt(2:5,:)'));
 
 
 k=kt;
 ppAnova=pp
 
 
[results, means, h, nms]=multcompare(stats);
tbl = array2table(results,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
clusters=statGroupClusters(results,means) 
   writetable(tbl,'AnovaMultCompLZDCRO.csv')

set(gca,'YScale','log')
title('LZD-CRO')
ylim([1.e-1 8])

end
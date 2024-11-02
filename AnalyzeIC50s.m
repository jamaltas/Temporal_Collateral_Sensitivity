function AnalyzeIC50s(csORcrflag,movMeanWindow,eps)
% look at the IC50 time series; csORcrflag=1 to plot prob of CR, 2 for prob
% of CS
%eps=0.3;  %threshold (remove data with abs < this)


edgevec4densityplot={linspace(-2,9,8),linspace(-6,6,10)};
%movMeanWindow is usually about 1-2



tmp=load('TempData.mat');
TempData=tmp.TempData;
TimeSeries=TempData.timeseries;
tvec=TempData.time;

% let's look at resistance to selecting drug vs collateral effects;  


resSel=NaN*TimeSeries(:,:,2:end);
resOther=TimeSeries(:,:,2:end);
    

for i=1:5  %5 drugs

    resOther([1:4]+(i-1)*4,i,:)=NaN;
    resSel([1:4]+(i-1)*4,i,:)=TimeSeries([1:4]+(i-1)*4,i,2:end);
  
        
end



resSel=resSel(isfinite(resSel));
resSel=reshape(resSel,20,1,[]);
resSel=repmat(resSel,1,5,1);

resOtherDiff=NaN*ones(size(TimeSeries));
resOtherDiff(:,:,2:4)=diff(resOther,1,3);
resOtherDiff(:,:,1)=resOther(:,:,1);
resOtherDiff=resOtherDiff(:,:,1:end-1);


resOtherBinary=NaN*zeros(size(resOther));
resOther(abs(resOther)<eps)=NaN;

if csORcrflag==1  %CR
resOtherBinary(resOther>0)=1;  %CR
resOtherBinary(resOther<=0)=-1;  %non-CR
plotbuffer=0.3;
ylab='Prob of CR';
ylimitProbPlot=[.2 1.1];
elseif csORcrflag==2  %CS
resOtherBinary(resOther<0)=1;  %CS
resOtherBinary(resOther>=0)=-1;  %non-CS
plotbuffer=0.55;
ylab='Prob of CS';
ylimitProbPlot=[0 .7];
end    
   
% % % figure
% % % for i=1:5
% % %  subplot(1,5,i)
% % % %     resStmp=resSel([1:4]+(i-1)*4,:,:);     %full values
% % % %     resOtmp=resOther([1:4]+(i-1)*4,:,:);
% % % %     
% % % 
% % %      resStmp=resSel([1:4]+(i-1)*4,:,:);     %binary
% % %      resOtmp=resOtherBinary([1:4]+(i-1)*4,:,:);
% % % 
% % %      [st,stind]=sort(resStmp(:)+1e-9*randn(size(resStmp(:))));
% % % 
% % %     plot(resStmp(:),resOtmp(:),'o')
% % %     hold on
% % %     plot(resStmp(stind),(1+smoothdata(resOtmp(stind),'movmean',movMeanWindow,'SamplePoints',st))/2,'-*k')
% % %     plot([0 8],[0 0],'--r')
% % %     
% % %     
% % % end


%figure
subplot(2,2,1)
title(strcat('\epsilon= ', ' ',num2str(eps,2)))
[st,stind]=sort(resSel(:)+1e-9*randn(size(resSel(:))));
hold on

%densityplot(resSel(:),resOther(:),'edges',edgevec4densityplot)
densityplot(resSel(:),resOther(:),'edges',edgevec4densityplot)
%scattercloud(resSel(:),resOther(:),25,2,'k.')
%densityplot(resSel(:),resOther(:))
colorbar off
hold on

%let's caluclate number of entries in each window
resSelsorted=resSel(stind);
resOthherBinarysorted=resOtherBinary(stind);
for i=1:length(resSelsorted)
   tmpct=find(resSelsorted<resSelsorted(i)+movMeanWindow & resSelsorted>resSelsorted(i)-movMeanWindow & isfinite(resSelsorted)& isfinite(resOthherBinarysorted));
    numNonzero(i)=numel(tmpct); 
end


scatter1 = scatter(resSel(:),resOther(:),5,'MarkerFaceColor','k','MarkerEdgeColor','k'); 
scatter1.MarkerFaceAlpha = .15;
scatter1.MarkerEdgeAlpha = .15;
%plot(resSel(:),resOther(:),'.k','markerfacecolor','k')
mean_ydata=(smoothdata(resOther(stind),'movmean',movMeanWindow,'SamplePoints',st));
std_ydata=movstd(resOther(stind),movMeanWindow,'omitnan','SamplePoints',st);
shadedErrorBar(resSelsorted,mean_ydata,std_ydata./sqrt(numNonzero'),'lineprops',{'k-','linewidth',2})
%plot(resSelsorted,mean_ydata,'-k','linewidth',2)

cm=flipud(gray(256));
colormap(cm)


plot([0 8],[0 0],'--r','linewidth',2)
axis square
xlabel('Resistance to Selecting Drug')
ylabel('Collateral Effects')
set(gca,'fontsize',12)
%ylim([-3 3])
xlim([0 7.5])
ylim([-1.5 2.5])


subplot(2,2,3)
[st,stind]=sort(resSel(:)+1e-9*randn(size(resSel(:))));
hold on

densityplot(resSel(:),resOtherDiff(:),'edges',edgevec4densityplot)
%densityplot(resSel(:),resOtherDiff(:))
colorbar off
hold on
scatter1 = scatter(resSel(:),resOtherDiff(:),5,'MarkerFaceColor','k','MarkerEdgeColor','k'); 
scatter1.MarkerFaceAlpha = .15;
scatter1.MarkerEdgeAlpha = .15;


mean_ydata=(smoothdata(resOtherDiff(stind),'movmean',movMeanWindow,'SamplePoints',st));
std_ydata=movstd(resOtherDiff(stind),movMeanWindow,'omitnan','SamplePoints',st);
shadedErrorBar(resSelsorted,mean_ydata,std_ydata./sqrt(numNonzero'))

plot(resSelsorted,mean_ydata,'-k','linewidth',2)

plot([0 8],[0 0],'--r','linewidth',2)
axis square
xlabel('Resistance to Selecting Drug')
ylabel('Inst. Collateral Effects')
set(gca,'fontsize',12)
ylim([-1.5 1.5])
xlim([0 7.5])




subplot(2,2,4)
d1tm=resSel(stind);
d2tm=resOtherBinary(stind);
% d1tm=d1tm(d2tm>0);
% d2tm=d2tm(d2tm>0);

plot(d1tm, d2tm/50+plotbuffer,'.k','markersize',4)
hold on
ylim(ylimitProbPlot)
xlim([0 7.1])
 d1tm=d1tm(d2tm>0);
 d2tm=d2tm(d2tm>0);
plot(d1tm, d2tm/50+plotbuffer,'.k','markersize',4)



[st,stind]=sort(resSel(:)+1e-9*randn(size(resSel(:))));
hold on
p1tmp=(1+smoothdata(resOtherBinary(stind),'movmean',movMeanWindow,'SamplePoints',st))/2;





stderr=movstd(resOtherBinary(stind),movMeanWindow,'omitnan','SamplePoints',st);
stderr=stderr./(2*sqrt(numNonzero'));
p2tmp=1-p1tmp;
%plot(resSel(stind),p1tmp,'ok','markersize',6,'markerfacecolor','k')

shadedErrorBar(resSelsorted,p1tmp,stderr)

axis square
set(gca,'fontsize',12)
xlabel('Resistance to Selecting Drug')
ylabel(ylab)

%[p1,h1]=corr(resSel(isfinite(resSel)&isfinite(resOther)),resOther(isfinite(resSel)&isfinite(resOther)))

resSelLRtmp=resSel(stind);
resOtherBinaryLRtmp=resOtherBinary(stind);

resSelLR=resSel(stind);
resOtherBinaryLR=resOtherBinary(stind);
resOtherBinaryLR(resOtherBinaryLR==-1)=2;  % make them positive integers



resSelLR=resSelLR(isfinite(resOtherBinaryLRtmp) & isfinite(resSelLRtmp));
resOtherBinaryLR=resOtherBinaryLR(isfinite(resOtherBinaryLRtmp) & isfinite(resSelLRtmp));
unique(resOtherBinaryLR)

[B,dev,stats] = mnrfit(resSelLR(:),resOtherBinaryLR(:))

statsIs=stats.p

LL = stats.beta - 1.96.*stats.se
UL = stats.beta + 1.96.*stats.se

subplot(2,2,2)
xx=linspace(0,8,100);
plot(resSel(stind),log(p1tmp./p2tmp),'ko')
hold on
plot(xx,B(1)+B(2)*xx,'r')
axis square
xlabel('Resistance to Selecting Drug')
ylabel('ln(p/1-p)')
set(gca,'fontsize',12)
xlim([0 8])


% subplot(2,2,4)
% plot(resSelsorted,numNonzero,'o')
% xlabel('Resistance to Selecting Drug')
% ylabel('Elements in window')
% axis square
function out=MakeHeterogeneityFigurePanels
% this makes a figure using Anh's new data on heterogeneity (in response to
% reviewer requests)


function cstmp=plotProb(outtmp)
ktmp=outtmp.k;
kerrtmp=outtmp.kerr;
wttmp=outtmp.wtMean;
wterrtmp=outtmp.wtStErr;

cstmp=zeros(size(ktmp));
% count as CS/CR if std errorbars don't overlap
% cstmp(ktmp+1*kerrtmp<wttmp-1*wterrtmp)=-1;
% cstmp(ktmp-1*kerrtmp>wttmp+1*wterrtmp)=+1;

% count as CS/CR if outside of 3 wt STD
cstmp(ktmp<wttmp-3*wterrtmp)=-1;
cstmp(ktmp>wttmp+3*wterrtmp)=+1;


cstmpcopy=cstmp;

cstmp=cstmp(2:end,:);  % don't need the WT;
Ntmp=size(cstmp,2);
probCR=sum(cstmp==1,2)/Ntmp;
probCRerr=sqrt(probCR.*(1-probCR)./Ntmp);
probCS=sum(cstmp==-1,2)/size(cstmp,2);
probCSerr=sqrt(probCS.*(1-probCS)./Ntmp);


probCRWT=sum(cstmpcopy==1,2)/Ntmp;
probCSWT=sum(cstmpcopy==-1,2)/Ntmp;

% CR_falsepos=probCRWT(1)
% CS_falsepos=probCSWT(1)



errorbar([2:2:8],probCR,probCRerr,'-ro','markerfacecolor','r','MarkerSize',7)
hold on
errorbar([2:2:8],probCS,probCSerr,'-ob','markerfacecolor','b','MarkerSize',7)
axis square
xlabel('Days')
ylabel('Freq of Effect')
xlim([1 9])
ylim([-0.05 1.05])
set(gca,'fontsize',12)
end

currd=pwd;


figure

cd LINE-CRO/expt_3/16hrs/
subplot(2,3,1)
[outLZDCRO]=QuickLook2LineCROFigPanel;
box on
cd(currd)


cd CRO-DOX/matFiles
subplot(2,3,2)
[outCRODOX]=QuickLook2CRODOXFigPanel;
box on
cd(currd)

cd CIP_CRO/16hrs
subplot(2,3,3)
[outCIPCRO]=QuickLook2CIPCroFigPanel;
box on
cd(currd)

out.outCIPCRO=outCIPCRO;
out.outCRODOX=outCRODOX;
out.outLZDCRO=outLZDCRO;


% now let's look for CS and CR that is statistically sign and plot probs
days=repmat([2:2:8]',1,12);

subplot(2,3,4)
disp('LZD-CRO Log Regs')
csMat=plotProb(outLZDCRO);
csMat2=csMat;
csMat(csMat==0 | csMat==1)=2;
csMat(csMat==-1)=1;
[B,dev,stats]=mnrfit(days(:),csMat(:));
betaCS=stats.beta
pvalCS=stats.p


csMat2(csMat2==0 | csMat2==-1)=2;
csMat2(csMat2==1)=1;
[B,dev,stats]=mnrfit(days(:),csMat2(:));
betaCR=stats.beta
pvalCR=stats.p
%%%%%%%%%%%%%
%%%%%%%%%%%%%

subplot(2,3,5)
disp('CRO-DOX Log Regs')
csMat=plotProb(outCRODOX);
csMat2=csMat;
csMat(csMat==0 | csMat==1)=2;
csMat(csMat==-1)=1;
[B,dev,stats]=mnrfit(days(:),csMat(:));
betaCS=stats.beta
pvalCS=stats.p


csMat2(csMat2==0 | csMat2==-1)=2;
csMat2(csMat2==1)=1;
[B,dev,stats]=mnrfit(days(:),csMat2(:));
betaCR=stats.beta
pvalCR=stats.p



subplot(2,3,6)
disp('CIP-CRO Log Regs')
csMat=plotProb(outCIPCRO);
csMat2=csMat;
csMat(csMat==0 | csMat==1)=2;
csMat(csMat==-1)=1;
[B,dev,stats]=mnrfit(days(:),csMat(:));
betaCS=stats.beta
pvalCS=stats.p


csMat2(csMat2==0 | csMat2==-1)=2;
csMat2(csMat2==1)=1;
[B,dev,stats]=mnrfit(days(:),csMat2(:));
betaCR=stats.beta
pvalCR=stats.p
%%%%%%%%%%%

disp('Comment out Anova Analysis in QuickLook2...FigPanel.m files to get pretty figures')


end

px=200;py=200;tlimit=240;patch=px*py;w=200;
col=find(sum(invasion_sto_A)==0);
invasion_sto_A(:,col)=[];
ww=w-length(col);
for y=1:ww
Map= invasion_sto_A(:,y);
invasionmap=zeros(px,py);
 for i = 1:px
   for j = 1:py
       invasionmap(i,j)= Map((j-1)*px+i);
   end
 end

cutoff_time=0:max(Map);
Invaded_patch=zeros(numel(cutoff_time),1);
%Invaded_fraction=zeros(numel(cutoff_time),1);
fraction_map=zeros(px,py);
for jj=1:numel(cutoff_time)
 for i = 1:px
   for j = 1:py
      if  invasionmap(i,j)<=cutoff_time(jj)
        fraction_map(i,j)=1;
      else fraction_map(i,j)=0;
      end
   end
   end
 Invaded_patch(jj)=sum(fraction_map,'all');
 Invaded_fraction(1,jj)=Invaded_patch(jj)/patch;

end
for i=1:min(cutoff_time)-1
    Invaded_fraction(1,i)=0;
end
for i=numel(cutoff_time)+1:tlimit
    Invaded_fraction(1,i)=1;
end
Invaded_fraction2(y,:)=Invaded_fraction;
clear Invaded_fraction
end

figure(20)
h1=plot([1:240],mean(Invaded_fraction2(:,1:240)),'LineWidth',2,'color','k');
hold on
A=sort(Invaded_fraction2(:,1:240));
A([1:round(ww*.025),round(ww*.975):ww],:)=[];
h2=plot([1:240],max(A),'LineWidth',2,'color', 'b');
hold on
h3=plot([1:240],min(A),'LineWidth',2,'color', 'r');
xlabel('Time (in months)','FontSize',14,'FontName','Times New Roman','FontWeight','bold','Color','k')
ylabel('Invaded fraction of the landscape by adults','FontSize',14,'FontName','Times New Roman','FontWeight','bold','Color','k')
legend([h1 h2 h3],{'Stochastic mean','Upper limit of 95% CI','Lower limit of 95% CI'});

xticks(0:30:240);
xlim([0 240])
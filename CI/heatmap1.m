tic
clc;
KK=3.4;%invasion threshold.
tlimit=240;%simulation time up to tlimit
px=20;py=20; %200*200 squared cell
patch=px*py;%No of pacthes
k=4*patch; %no of stages
ddL=0;% dispersal distance for L and N due to tick movement
ddA=0;% dispersal distance for A due to tick movement
ph_L=1;% dispersal distance for L due to host movement
ph_A=22;% dispersal distance for A due to host movement
nzero=zeros(4*patch,1)';
sx=1;sy=1; %initial position(left top)
StartingPatch=(sy-1)*px+sx; 
nzero(StartingPatch)=0;
nzero(StartingPatch+patch)=0;
nzero(StartingPatch+2*patch)=0;
nzero(StartingPatch+3*patch)=6; %initial condition for all E,L,N,A

a=0;c=8.65; %Carrying capacity is 150;
h_L=1*5*3; %(fraction of the patches in the home range that a host visits in one unit of time*
%h_N=h_L; %number of patches in the home range*density for one patch)
h_A=1*500*0.0428; %we may want to do some sensitivity analysis on it and decide to change these values based on what we find.
y=0;
h_L=h_L*(1-y);
h_A=h_A*(1-y);

S_E=0.6198; %survival probabilities are same for all 12 months
S_L=0.5368;
S_N=0.5698;
S_A=0.7349;

gammaE=0.5952;
gammaL=(1).*[0 0 ones(1,6) zeros(1,4)]; %\gamma are not same
gammaN=(0.7881).*[0 0 0 ones(1,6) zeros(1,3)];
gammaA=(0.5630).*[ones(1,6) zeros(1,6)];


sigmaL=(1-exp(-0.1*h_L)).*[0 0 ones(1,6) zeros(1,4)];
sigmaN=(1-exp(-0.1*h_L)).*[0 0 0 ones(1,6) zeros(1,3)];
sigmaA=(1-exp(-2*h_A)).*[ones(1,6) zeros(1,6)];

beta=(2380.564).*[ones(1,6) zeros(1,6)];
conversion=1;






       TE=S_E*(1-gammaE);
       TE1=S_E*(gammaE);
       TL1=conversion*S_L.*gammaL.*sigmaL;
       TN1=conversion*S_N.*gammaN.*sigmaN;
       dd=conversion.*gammaA.*sigmaA;
       TL=S_L*(1-gammaL.*sigmaL);
       TN=S_N*(1-gammaN.*sigmaN);
       TA=S_A*(1-gammaA.*sigmaA);
       

   m_L15=distance_mat(2,ddL,px,py); %distance_mat function has to be in the directory
   m_A15=distance_mat(2,ddA,px,py);
   psi_L15=distance_mat(2,ph_L,px,py);
   psi_A15=distance_mat(2,ph_A,px,py);
  
  
invasion_temp=ones(patch,12);
invasion_det=zeros(patch,1,'single');
InvasionYes1=ones(patch,1,'single');
aaa=zeros(4*patch,tlimit,'single');
aaa(:,1)=nzero;
taaa1=zeros(patch,tlimit);
ttt=zeros(patch,tlimit);
%taaa1=zeros(12,patch);
 
  ff1=a.*gammaA.*sigmaA;
  ff2=c.*gammaA.*sigmaA;
  ff3=gammaA.*sigmaA;
  
  
  

  
         
 
   


w=10;%change it according to your needs. We used 200
invasion_sto=zeros(patch,w,'single');
tn=zeros(patch,w,tlimit,'single');
tt=zeros(patch,tlimit);

for y=1:w
invasion_temp2=ones(patch,12);
    fprintf('Current simulation number is %d\n',y)
    n= nzero;
  
    nprime=zeros(k,1,'single');
    nout=zeros(4*patch,tlimit,'single');
    InvasionYes=ones(patch,1,'single');
   tn(1,y,1)=nzero(3*patch+1);
   for t = 2:tlimit
  

   %fprintf('Current time number is %d\n',t)
    jj=mod(t,12);
    if jj==0
        jj=12;
    end
   
   m_L15=distance_mat(2,ddL,px,py);
   m_N15=m_L15;
   m_A15=distance_mat(2,ddA,px,py);
   
   psi_L15=distance_mat(2,ph_L,px,py);
   psi_N15=psi_L15;
   psi_A15=distance_mat(2,ph_A,px,py); 
   
   m_L15(:,3)=m_L15(:,3).*TL(jj);
   psi_L15(:,3)=psi_L15(:,3).*TL1(jj);
  
   psi_L15(:,2)=psi_L15(:,2)+patch;
   comb1=sortrows([m_L15;psi_L15]);
   mask2=mask_mat(comb1);
   comb12=comb1(:,2);
   comb13=comb1(:,3);
   m_N15(:,3)=m_N15(:,3).*TN(jj);
   psi_N15(:,3)=psi_N15(:,3).*TN1(jj);
   
   psi_N15(:,2)=psi_N15(:,2)+patch;
   comb2=sortrows([m_N15;psi_N15]);
    
   mask3=mask_mat(comb2);
   comb22=comb2(:,2);
   comb23=comb2(:,3);
   m_A15(:,3)=m_A15(:,3).*TA(jj);
   m_A15=sortrows(m_A15);
   mask4=mask_mat(m_A15);
   m_A152=m_A15(:,2);
   m_A153=m_A15(:,3);
   psi_A15(:,3)=psi_A15(:,3).*dd(jj);
   psi_A15=sortrows(psi_A15);
   psi_A15=psi_A15(:,3);
  

   
ff1=a.*gammaA.*sigmaA;
ff2=c.*gammaA.*sigmaA;
ff3=gammaA.*sigmaA;


rci=zeros(length(psi_A15),2,'single'); 
counter2=1;
for i = 1:px
    
   for j = 1:py
       
       for ss=max(1,i-ph_A):min(px,i+ph_A)
           for s=max(1,j-ph_A):min(py,j+ph_A)
               tmp=sqrt((i-ss)^2+(s-j)^2);
               if tmp<=ph_A
                  col=(j-1)*px+i;
                  row=(s-1)*px+ss;
                  rci(counter2,1)=row;
                  rci(counter2,2)=col;
                  counter2=counter2+1;
               end 
           end
       end  
   end
end
 rci=sortrows(rci);
 mask5=mask_mat(rci(:,1));
 rci2=rci(:,2);


  tmp8=[patch;2*patch;4*patch+1];
  tmp9=[patch+comb1(mask2(patch):length(comb1),2);4*patch+1];
  tmp10=[2*patch+comb2(mask3(patch):length(comb2),2);4*patch+1];
  tmp11=[3*patch+m_A15(mask4(patch):length(m_A15),2);4*patch+1];
  tmpc=3*patch+rci(:,1);
  tmp13=[rci2(mask5(patch):length(rci2))];
 

A=[TE;TE1;1-TE-TE1];

B=1-TL(jj)-TL1(jj);
C=1-TN(jj)-TN1(jj);
D=1-TA(jj);
x1=[comb13(mask2(patch):length(comb1));B]';
x2=[comb23(mask3(patch):length(comb2));C]';
x3=[m_A153(mask4(patch):length(m_A15));D]';
x5=(mask5(patch):length(rci));


      tmp1=zeros(3,patch);
      for i = 1:patch
      tmp1(:,i)=[i;patch+i;4*patch+1];
      end
      
      
tmp2=cell(1,patch);
tmp3=cell(1,patch);
tmp4=cell(1,patch);
tmp5=cell(1,patch);
tmp6=cell(1,patch);
tmp7=cell(1,patch);
tmp12=cell(1,patch);
for i=1:patch-1
tmp2(:,i)={[patch+comb12(mask2(i):mask2(i+1)-1);4*patch+1]};
tmp3(:,i)={[comb13(mask2(i):mask2(i+1)-1);B]'};
tmp4(:,i)={[2*patch+comb22(mask3(i):mask3(i+1)-1);4*patch+1]};
tmp5(:,i)={[comb23(mask3(i):mask3(i+1)-1);C]'};
tmp6(:,i)={[3*patch+m_A152(mask4(i):mask4(i+1)-1);4*patch+1]};
tmp7(:,i)={[m_A153(mask4(i):mask4(i+1)-1);D]'};
tmp12(:,i)={[rci2(mask5(i):mask5(i+1)-1)]};
end 
    
    
    
    
    
 
    
   nout(:,t)=n;
   trans1=zeros(k+1,1);
   trans2=zeros(k+1,1);
   trans3=zeros(k+1,1);
   trans4=zeros(k+1,1);


   tmpv=nout(tmpc,t);
   F3=(((beta.*(1+ff1(jj).*tmpv))./(1+ff2(jj).*tmpv+a.*(ff3(jj).*tmpv).^2)));
   F3=F3.*(psi_A15);

  

 for i = 1:patch-1
       if nout(i,t)>0
        trans1(tmp1(:,i))=trans1(tmp1(:,i))+mnrnd (nout(i,t),A')';
       end 
       if nout(patch+i,t)>0
        trans2(tmp2{i})=trans2(tmp2{i})+mnrnd(nout(patch+i,t),tmp3{i})';
       end
     if nout(2*patch+i,t)>0
        trans3(tmp4{i})=trans3(tmp4{i})+mnrnd (nout(2*patch+i,t),tmp5{i})';
     end
     if nout(3*patch+i,t)>0
       trans4(tmp6{i})=trans4(tmp6{i})+mnrnd(nout(3*patch+i,t),tmp7{i})';
     end
 end
  
 
    trans1(tmp8)=trans1(tmp8)+mnrnd (nout(patch,t),A')';
    trans2(tmp9)=trans2(tmp9)+mnrnd (nout(2*patch,t),x1)';
    trans3(tmp10)=trans3(tmp10)+mnrnd (nout(3*patch,t),x2)';   
    trans4(tmp11)=trans4(tmp11)+mnrnd(nout(4*patch,t),x3)';
    nprime=trans1+trans2+trans3+trans4;
    nprime=nprime(1:k);
    
 for i=1:patch-1
     if (nout(3*patch+i,t)>0)
nprime(tmp12{i})=nprime(tmp12{i})+poissrnd(F3(mask5(i):mask5(i+1)-1)*nout(3*patch+i,t))';
     end
 end 

   nprime(tmp13) = nprime(tmp13)+poissrnd(F3(x5)*nout(4*patch,t))';
   n = nprime;
   
      for i = 1:patch
            tn(i, y, t) = nprime((3 * patch + i));
        end

        % Invasion detection (after at least 12 time steps)
        if t >= 12
            for i = 1:patch
                if invasion_sto(i, y) == 0
                    start_idx = t - 11;
                    end_idx = t;
                    window_mean1 = mean(tn(i, y, start_idx:end_idx));
                    if window_mean1 >= 3.4
                        invasion_sto(i, y) = t;  % Store the time of invasion
                    end
                end
            end

            % Stop early if all patches have been invaded
            if all(invasion_sto(:, y) > 0)
                %fprintf('All patches invaded in simulation %d by time %d. Stopping early.\n', y, t);
                break;
            end
        end
    end
end


Map=round(mean(invasion_sto'));
invasionmap_sto=zeros(px,py);
 for i = 1:px
   for j = 1:py
       invasionmap_sto(i,j)= Map((j-1)*px+i);
   end
 end
 col=find(sum(invasion_sto)==0);
invasion_sto(:,col)=[];
ww=w-length(col);
for y=1:ww
Map= invasion_sto(:,y);
invasionmap=zeros(px,py);
 for i = 1:px
   for j = 1:py
       invasionmap(i,j)= Map((j-1)*px+i);
   end
 end

cutoff_time=0:max(max(invasionmap));
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
%legend([h1 h2 h3],{'Stochastic mean','Upper limit of 95% CI','Lower limit of 95% CI'});

xticks(0:30:240);
xlim([0 240])


toc
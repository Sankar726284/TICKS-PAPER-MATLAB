tic
clc;
KK=3.4;
tlimit=240;%simulation time up to tlimit
px=20;py=20;
patch=px*py;
k=4*patch;

ddL=0;% dispersal distance for L and N due to tick movement
ddA=0;% dispersal distance for A due to tick movement
ph_L=1;% dispersal distance for L due to host movement
ph_A=22;% dispersal distance for A due to host movement

nzero=zeros(4*patch,1)';

sx=1;sy=1;
StartingPatch=(sy-1)*px+sx; 
nzero(StartingPatch)=0;
nzero(StartingPatch+patch)=0;
nzero(StartingPatch+2*patch)=0;
nzero(StartingPatch+3*patch)=6;

a=0;c=8.65;

h_L=1*5*3;
h_A=1*500*0.042;
y=0;
h_L=h_L*(1-y);
h_A=h_A*(1-y);

S_E=0.6198;
S_L=0.5368;
S_N=0.5698;
S_A=0.7349;

gammaE=0.5952;
gammaL=1.*[0 0 ones(1,6) zeros(1,4)];
sigmaL=(1-exp(-0.1*h_L)).*[0 0 ones(1,6) zeros(1,4)];

gammaN=0.7881.*[0 0 0 ones(1,6) zeros(1,3)];
sigmaN=(1-exp(-0.1*h_L)).*[0 0 0 ones(1,6) zeros(1,3)];
gammaA=0.5630.*[ones(1,6) zeros(1,6)];
sigmaA=(1-exp(-2.0*h_A)).*[ones(1,6) zeros(1,6)];
beta=(2380.564).*[ones(1,6) zeros(1,6)];
conversion=1;











S_E2=S_E*0.5; %50 percent reduction in survival
S_L2=S_L*0.5;
S_N2=S_N*0.5;
S_A2=S_A*0.5;

S_E2=S_E*0.25;  %75 percent reduction in survival
S_L2=S_L*0.25;
S_N2=S_N*0.25;
S_A2=S_A*0.25;



       TE=S_E*(1-gammaE);
       TE1=S_E*(gammaE);
       TL1=conversion*S_L.*gammaL.*sigmaL;
       TN1=conversion*S_N.*gammaN.*sigmaN;
       dd=conversion.*gammaA.*sigmaA;
       TL=S_L*(1-gammaL.*sigmaL);
       TN=S_N*(1-gammaN.*sigmaN);
       TA=S_A*(1-gammaA.*sigmaA);
       
       TE2=S_E2*(1-gammaE);
       TE12=S_E2*(gammaE);
       TL12=conversion*S_L2.*gammaL.*sigmaL;
       TN12=conversion*S_N2.*gammaN.*sigmaN;
       TL2=S_L2*(1-gammaL.*sigmaL);
       TN2=S_N2*(1-gammaN.*sigmaN);
       TA2=S_A2*(1-gammaA.*sigmaA);
       
     dd1=zeros(patch);
     gg=zeros(1,patch);

 for i = 1:px
   for j = 1:py
      
               if 2<sqrt((1-i)^2+(1-j)^2) && sqrt((1-i)^2+(1-j)^2)<=7
                 
                  dd1(i,j)= 1;
                  gg((j-1)*px+i)=(j-1)*px+i;
               end 
  end
end   
    nn=find(gg>0);   
 
   m_L15=distance_mat(2,ddL,px,py);
   m_A15=distance_mat(2,ddA,px,py);
   psi_L15=distance_mat(2,ph_L,px,py);
   psi_A15=distance_mat(2,ph_A,px,py);

   BB=[];
   CC=[];
   for i=1:length(nn)
   m_L_rows1(:,i)=find((m_L15(:,1)==nn(i)));
   %m_L_rows2(:,i)=find((m_L15(:,1)~=nn(i)));
   end
   
   for i=1:length(nn)
   psi_L_rows1=find((psi_L15(:,1)==nn(i)));
   BB=[BB;psi_L_rows1];
   
   end
  
 
   for i=4:14
      m_L15(:,i)=m_L15(:,3);
      psi_L15(:,i)=psi_L15(:,3);
      m_A15(:,i)=m_A15(:,3);
      psi_A15(:,i)=psi_A15(:,3);
   end
   
   m_N15=m_L15;
   psi_N15=psi_L15;
   
   
for jj=1:12 %12
    for i=1:length(m_L15)
      if any(i==m_L_rows1)
    m_L15([i],jj+2)=m_L15([i],jj+2)*TL2(jj);
    m_N15([i],jj+2)=m_N15([i],jj+2).*TN2(jj);
    m_A15([i],jj+2)=m_A15([i],jj+2)*TA2(jj);
      else
     m_L15([i],jj+2)=m_L15([i],jj+2)*TL(jj);
     m_N15([i],jj+2)=m_N15([i],jj+2).*TN(jj);
     m_A15([i],jj+2)=m_A15([i],jj+2)*TA(jj);
      end
    end
    
    for i=1:length(psi_L15)
      if any(i==BB)
    psi_L15([i],jj+2)=psi_L15([i],jj+2)*TL12(jj);
    psi_N15([i],jj+2)=psi_N15([i],jj+2).*TN12(jj);
    
      else
    psi_L15([i],jj+2)=psi_L15([i],jj+2)*TL1(jj);
    psi_N15([i],jj+2)=psi_N15([i],jj+2).*TN1(jj);
      end
    end

    
    
    
     
    psi_A15(:,jj+2)=psi_A15(:,jj+2)*dd(jj);
end
 
  %psi_L15
  
   psi_L15(:,2)=psi_L15(:,2)+patch;
   comb1=sortrows([m_L15;psi_L15]);
   
   mask2=mask_mat(comb1);
   psi_N15(:,2)=psi_N15(:,2)+patch;
   comb2=sortrows([m_N15;psi_N15]);
   mask3=mask_mat(comb2);
   m_A15=sortrows(m_A15);
   mask4=mask_mat(m_A15);   
   psi_A15=sortrows(psi_A15);
  
  
   
 

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
  tmp11=[3*patch+m_A15(mask4(patch):length(m_A15(:,1)),2);4*patch+1];
  tmpc=3*patch+rci(:,1);
  tmp13=[rci2(mask5(patch):length(rci2))];




A1=[TE;TE1;1-TE-TE1];
A2=[TE2;TE12;1-TE2-TE12];

     tmp1=zeros(3,patch);
      for i = 1:patch
      tmp1(:,i)=[i;patch+i;4*patch+1];
      end
      
      
tmp2=cell(patch,1);
tmp3=cell(patch,1,12);
tmp4=cell(patch,1);
tmp5=cell(patch,1,12);
tmp6=cell(patch,1);
tmp7=cell(patch,1,12);
tmp12=cell(patch,1);
x1=cell(12,1);
x2=cell(12,1);
x3=cell(12,1);

for jj=1:12
  K1=comb1(mask2(patch):length(comb1(:,1)),jj+2);
  K2=comb2(mask3(patch):length(comb2(:,1)),jj+2);
  K3=m_A15(mask4(patch):length(m_A15(:,1)),jj+2);
x1(jj,:)={[K1;1-sum(K1)]'};
x2(jj,:)={[K2;1-sum(K2)]'};
x3(jj,:)={[K3;1-sum(K3)]'};
  
for i=1:patch-1
tmp2(i,:)={[patch+comb1(mask2(i):mask2(i+1)-1,2);4*patch+1]};
K4=comb1(mask2(i):mask2(i+1)-1,jj+2);
tmp3(i,:,jj)={[K4;1-sum(K4)]'};
tmp4(i,:)={[2*patch+comb2(mask3(i):mask3(i+1)-1,2);4*patch+1]};
K5=comb2(mask3(i):mask3(i+1)-1,jj+2);
tmp5(i,:,jj)={[K5;1-sum(K5)]'};
tmp6(i,:)={[3*patch+m_A15(mask4(i):mask4(i+1)-1,2);4*patch+1]};
K6=m_A15(mask4(i):mask4(i+1)-1,jj+2);
tmp7(i,:,jj)={[K6;1-sum(K6)]'};
tmp12(i,:)={[rci2(mask5(i):mask5(i+1)-1)]};
end 
end

x5=(mask5(patch):length(rci));


w=5;
Invaded_fraction=zeros(1,tlimit);
Invaded_fraction2=zeros(w,tlimit);
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
   for t = 2:tlimit-1
  

   %fprintf('Current time number is %d\n',t)
    jj=mod(t,12);
    if jj==0
        jj=12;
    end
   
   nout(:,t)=n;
   trans1=zeros(k+1,1);
   trans2=zeros(k+1,1);
   trans3=zeros(k+1,1);
   trans4=zeros(k+1,1);


   tmpv=nout(tmpc,t);
   F3=(((beta.*(1+ff1(jj).*tmpv))./(1+ff2(jj).*tmpv+a.*(ff3(jj).*tmpv).^2)));
   F3=F3.*(psi_A15(:,jj+2));
   
   for i = 1:patch-1
       if nout(i,t)>0 
           if any(i~=nn)
               tmp1=[i;patch+i;4*patch+1];
               trans1(tmp1)=trans1(tmp1)+mnrnd (nout(i,t),A1')';
           else
               tmp1=[i;patch+i;4*patch+1];
               trans1(tmp1)=trans1(tmp1)+mnrnd (nout(i,t),A2')';
           end
       end
       if nout(patch+i,t)>0
        trans2(tmp2{i,:})=trans2(tmp2{i,:})+mnrnd(nout(patch+i,t),tmp3{i,:,jj})';
       end
     if nout(2*patch+i,t)>0
        trans3(tmp4{i,:})=trans3(tmp4{i,:})+mnrnd (nout(2*patch+i,t),tmp5{i,:,jj})';
     end
 
     if nout(3*patch+i,t)>0
       trans4(tmp6{i,:})=trans4(tmp6{i,:})+mnrnd(nout(3*patch+i,t),tmp7{i,:,jj})';
     end
   end
 
    trans1(tmp8)=trans1(tmp8)+mnrnd (nout(patch,t),A1')';
    trans2(tmp9)=trans2(tmp9)+mnrnd (nout(2*patch,t),x1{jj,:})';
    trans3(tmp10)=trans3(tmp10)+mnrnd (nout(3*patch,t),x2{jj,:})';   
    trans4(tmp11)=trans4(tmp11)+mnrnd(nout(4*patch,t),x3{jj,:})';
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

h=zeros(w,240);
for y=1:w
Map= invasion_sto(:,y);
invasionmap=zeros(px,py);
 for i = 1:px
   for j = 1:py
       invasionmap(i,j)= Map((j-1)*px+i);
   end
 end
for t=1:240
h(y,t)=length(find(invasionmap>0 & invasionmap<=t))/patch;
end
end
plot(mean(h),'LineWidth',3)
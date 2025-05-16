tic
clc;
KK=3.4;
tlimit=240;%simulation time up to tlimit
ll=2;
total_det=zeros(1,ll+1);
px=1;py=1;
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

a=0;c=8.65; %Carrying capacity is 150;  146(1%), 151(10%)

h_L=1*5*3; %(fraction of the patches in the home range that a host visits in one unit of time*
%h_N=h_L; %number of patches in the home range*density for one patch)
h_A=1*500*0.042; %we may want to do some sensitivity analysis on it and decide to change these values based on what we find.
y=0;
h_L=h_L*(1-y);
h_A=h_A*(1-y);

S_E=0.6198;
gammaE=0.5952;

S_L=0.5368;
gammaL=1.*[0 0 ones(1,6) zeros(1,4)];
sigmaL=(1-exp(-0.1*h_L)).*[0 0 ones(1,6) zeros(1,4)];

S_N=0.5698;
gammaN=0.7881.*[0 0 0 ones(1,6) zeros(1,3)];
sigmaN=(1-exp(-0.1*h_L)).*[0 0 0 ones(1,6) zeros(1,3)];

S_A=0.7349;
gammaA=0.5630.*[ones(1,6) zeros(1,6)];
sigmaA=(1-exp(-2.0*h_A)).*[ones(1,6) zeros(1,6)];
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
       

   m_L15=distance_mat(2,ddL,px,py);
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
  
  
  
for t=1:tlimit
 fprintf('Current simulation(T) number is %d\n',t)  
    jj=mod(t,12);
    if jj==0
        jj=12;
    end
  
ff=zeros(length(psi_A15(:,1)),3,'single');
counter2=1;

for i = 1:px
    
   for j = 1:py
       
       for ss=max(1,i-ph_A):min(px,i+ph_A)
           for s=max(1,j-ph_A):min(py,j+ph_A)
               tmp=sqrt((i-ss)^2+(s-j)^2);
               if tmp<=ph_A
                  row=(s-1)*px+ss;
                  col=(j-1)*px+i;
                  ff(counter2,1)=row;
                  ff(counter2,2)=col;
                  
                  ff(counter2,3)=(((beta(jj)*(1+ff1(jj)*aaa(3*patch+col,t)))/(1+ff2(jj)*aaa(3*patch+col,t)+a*(ff3(jj)*aaa(3*patch+col,t))^2)));
                  
                 counter2=counter2+1;
               end 
           end
       end  
   end
end

 ff=sortrows(ff);
 ff=ff(:,3).*(dd(jj).*psi_A15(:,3));
 
  
   
   mask1=mask_mat(psi_A15);
   mask2=mask_mat(m_L15);
   mask3=mask_mat(psi_L15);
   mask4=mask_mat(m_A15);
  
  
         
 
     for i=1:patch-1
         aaa(i,t+1)=TE*aaa(i,t)+ ff(mask1(i):mask1(i+1)-1,1)'*aaa((3*patch+[psi_A15(mask1(i):mask1(i+1)-1,2)]),t);
         aaa(patch+i,t+1)=TE1*aaa(i,t)+TL(jj)*m_L15(mask2(i):mask2(i+1)-1,3)'*aaa((patch+[m_L15(mask2(i):mask2(i+1)-1,2)]),t);
         aaa(2*patch+i,t+1)=TL1(jj)*psi_L15(mask3(i):mask3(i+1)-1,3)'*aaa((patch+[psi_L15(mask3(i):mask3(i+1)-1,2)]),t)+TN(jj)*m_L15(mask2(i):mask2(i+1)-1,3)'*aaa((2*patch+[m_L15(mask2(i):mask2(i+1)-1,2)]),t);
         aaa(3*patch+i,t+1)=TN1(jj)*psi_L15(mask3(i):mask3(i+1)-1,3)'*aaa((2*patch+[psi_L15(mask3(i):mask3(i+1)-1,2)]),t)+TA(jj)*m_A15(mask4(i):mask4(i+1)-1,3)'*aaa((3*patch+[m_A15(mask4(i):mask4(i+1)-1,2)]),t);       
     end
              
         aaa(patch,t+1)=TE*aaa(patch,t)+ff(mask1(patch):length(psi_A15(:,1)),1)'*aaa((3*patch+[psi_A15(mask1(patch):length(psi_A15(:,1)),2)]),t);
         aaa(patch+patch,t+1)=TE1*aaa(patch,t)+TL(jj)*m_L15(mask2(patch):length(m_L15(:,1)),3)'*aaa((patch+[m_L15(mask2(patch):length(m_L15(:,1)),2)]),t);
         aaa(2*patch+patch,t+1)=TL1(jj)*psi_L15(mask3(patch):length(psi_L15(:,1)),3)'*aaa((patch+[psi_L15(mask3(patch):length(psi_L15(:,1)),2)]),t)+TN(jj)*m_L15(mask2(patch):length(m_L15(:,1)),3)'*aaa((2*patch+[m_L15(mask2(patch):length(m_L15(:,1)),2)]),t);
         aaa(3*patch+patch,t+1)=TN1(jj)*psi_L15(mask3(patch):length(psi_L15(:,1)),3)'*aaa((2*patch+[psi_L15(mask3(patch):length(psi_L15(:,1)),2)]),t)+TA(jj)*m_A15(mask4(patch):length(m_A15(:,1)),3)'*aaa((3*patch+[m_A15(mask4(patch):length(m_A15(:,1)),2)]),t);
        
       
     
   for i = 1:patch
  taaa1(i, t) = aaa(3 * patch + i, t);
  end
    
  invasion_det = zeros(patch,1); % To store the value of x for each row

for i = 1:patch  % Loop through each row
    start_idx = 1;
    end_idx = 12;
    
    while end_idx <= tlimit  % Ensure the window stays within bounds
        window_mean = mean(taaa1(i, start_idx:end_idx)); % Compute mean
        
        if window_mean >= 3.4
            invasion_det(i) = end_idx; % Store the index where condition is met
            break;
        else
            % Shift the window by removing the first element and adding the next
            start_idx = start_idx + 1;
            end_idx = end_idx + 1;
        end
    end
end
if all(invasion_det > 0)
    fprintf('Invasion detected in all patches. Stopping simulation early at t = %d\n', t);
    break;
end
end



invasionmap_det=zeros(px,py,'single');
 for i = 1:px
   for j = 1:py
       invasionmap_det(i,j)= invasion_det((j-1)*px+i);
   end
 end
        
total_det(1)=max(max(invasionmap_det));

for l=1:ll
px=20*l;py=20*l;
patch=px*py;
k=4*patch;
nzero=zeros(4*patch,1)';
sx=1;sy=1;
StartingPatch=(sy-1)*px+sx; 
nzero(StartingPatch)=0;
nzero(StartingPatch+patch)=0;
nzero(StartingPatch+2*patch)=0;
nzero(StartingPatch+3*patch)=6;

       

   m_L15=distance_mat(2,ddL,px,py);
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
 
  ff1=a.*gammaA.*sigmaA;
  ff2=c.*gammaA.*sigmaA;
  ff3=gammaA.*sigmaA;
  
  
  
for t=1:tlimit
 %fprintf('Current simulation(T) number is %d\n',t)  
    jj=mod(t,12);
    if jj==0
        jj=12;
    end
  
ff=zeros(length(psi_A15(:,1)),3,'single');
counter2=1;

for i = 1:px
    
   for j = 1:py
       
       for ss=max(1,i-ph_A):min(px,i+ph_A)
           for s=max(1,j-ph_A):min(py,j+ph_A)
               tmp=sqrt((i-ss)^2+(s-j)^2);
               if tmp<=ph_A
                  row=(s-1)*px+ss;
                  col=(j-1)*px+i;
                  ff(counter2,1)=row;
                  ff(counter2,2)=col;
                  
                  ff(counter2,3)=(((beta(jj)*(1+ff1(jj)*aaa(3*patch+col,t)))/(1+ff2(jj)*aaa(3*patch+col,t)+a*(ff3(jj)*aaa(3*patch+col,t))^2)));
                  
                 counter2=counter2+1;
               end 
           end
       end  
   end
end

 ff=sortrows(ff);
 ff=ff(:,3).*(dd(jj).*psi_A15(:,3));
 
  
   
   mask1=mask_mat(psi_A15);
   mask2=mask_mat(m_L15);
   mask3=mask_mat(psi_L15);
   mask4=mask_mat(m_A15);
  
     for i=1:patch-1
         aaa(i,t+1)=TE*aaa(i,t)+ ff(mask1(i):mask1(i+1)-1,1)'*aaa((3*patch+[psi_A15(mask1(i):mask1(i+1)-1,2)]),t);
         aaa(patch+i,t+1)=TE1*aaa(i,t)+TL(jj)*m_L15(mask2(i):mask2(i+1)-1,3)'*aaa((patch+[m_L15(mask2(i):mask2(i+1)-1,2)]),t);
         aaa(2*patch+i,t+1)=TL1(jj)*psi_L15(mask3(i):mask3(i+1)-1,3)'*aaa((patch+[psi_L15(mask3(i):mask3(i+1)-1,2)]),t)+TN(jj)*m_L15(mask2(i):mask2(i+1)-1,3)'*aaa((2*patch+[m_L15(mask2(i):mask2(i+1)-1,2)]),t);
         aaa(3*patch+i,t+1)=TN1(jj)*psi_L15(mask3(i):mask3(i+1)-1,3)'*aaa((2*patch+[psi_L15(mask3(i):mask3(i+1)-1,2)]),t)+TA(jj)*m_A15(mask4(i):mask4(i+1)-1,3)'*aaa((3*patch+[m_A15(mask4(i):mask4(i+1)-1,2)]),t);       
     end
              
         aaa(patch,t+1)=TE*aaa(patch,t)+ff(mask1(patch):length(psi_A15(:,1)),1)'*aaa((3*patch+[psi_A15(mask1(patch):length(psi_A15(:,1)),2)]),t);
         aaa(patch+patch,t+1)=TE1*aaa(patch,t)+TL(jj)*m_L15(mask2(patch):length(m_L15(:,1)),3)'*aaa((patch+[m_L15(mask2(patch):length(m_L15(:,1)),2)]),t);
         aaa(2*patch+patch,t+1)=TL1(jj)*psi_L15(mask3(patch):length(psi_L15(:,1)),3)'*aaa((patch+[psi_L15(mask3(patch):length(psi_L15(:,1)),2)]),t)+TN(jj)*m_L15(mask2(patch):length(m_L15(:,1)),3)'*aaa((2*patch+[m_L15(mask2(patch):length(m_L15(:,1)),2)]),t);
         aaa(3*patch+patch,t+1)=TN1(jj)*psi_L15(mask3(patch):length(psi_L15(:,1)),3)'*aaa((2*patch+[psi_L15(mask3(patch):length(psi_L15(:,1)),2)]),t)+TA(jj)*m_A15(mask4(patch):length(m_A15(:,1)),3)'*aaa((3*patch+[m_A15(mask4(patch):length(m_A15(:,1)),2)]),t);
        
       
     
     
   for i = 1:patch
    taaa1(i, t) = aaa(3 * patch + i, t);
   end
    
invasion_det = zeros(patch,1); % To store the value of x for each row

for i = 1:patch  % Loop through each row
    start_idx = 1;
    end_idx = 12;
    
    while end_idx <= tlimit  % Ensure the window stays within bounds
        window_mean = mean(taaa1(i, start_idx:end_idx)); % Compute mean
        
        if window_mean >= 3.4
            invasion_det(i) = end_idx; % Store the index where condition is met
            break;
        else
            % Shift the window by removing the first element and adding the next
            start_idx = start_idx + 1;
            end_idx = end_idx + 1;
        end
    end
end
if all(invasion_det > 0)
    fprintf('Invasion detected in all patches. Stopping simulation early at t = %d\n', t);
    break;
end
end

invasionmap_det=zeros(px,py,'single');
 for i = 1:px
   for j = 1:py
       invasionmap_det(i,j)= invasion_det((j-1)*px+i);
   end
 end
 
 total_det(l+1)=max(max(invasionmap_det));
end
total_det








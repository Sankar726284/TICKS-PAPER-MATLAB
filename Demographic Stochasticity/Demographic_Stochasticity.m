tic
clc;
KK=3.4;
tlimit=240;%simulation time up to tlimit
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
a=0;c=8.65; %Carrying capacity is 150;
%a=0;c=1.29; %Carrying capacity is 1000;



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

dd1=zeros(patch);
dd2=zeros(patch);
dd3=zeros(patch);
for i = 1:px
   for j = 1:py
       for ss=max(1,i-ddL):min(px,i+ddL)
           for s=max(1,j-ddL):min(py,j+ddL)
               if sqrt((i-ss)^2+(s-j)^2)<=ddL
                  deltaL=2;%trucated discrete Laplace with a scale parameter deltaL
                  dd1((s-1)*px+ss,(j-1)*px+i)= 1/(2*deltaL)*exp(-sqrt((i-ss)^2+(s-j)^2)/deltaL);
               end 
           end
       end  
       
       for ss=max(1,i-ddA):min(px,i+ddA)
           for s=max(1,j-ddA):min(py,j+ddA)
               if sqrt((i-ss)^2+(s-j)^2)<=ddA
                  deltaA=2;%trucated discrete Laplace with a scale parameter deltaA
                  dd3((s-1)*px+ss,(j-1)*px+i)= 1/(2*deltaA)*exp(-sqrt((i-ss)^2+(s-j)^2)/deltaA);
               end    
           end
       end 
   end
end
m_L=(dd1./sum(dd1));
m_N=m_L;
m_A=(dd3./sum(dd3));

psi1=zeros(patch);
psi2=zeros(patch);
for i = 1:px
   for j = 1:py
       for ss=max(1,i-ph_L):min(px,i+ph_L)
           for s=max(1,j-ph_L):min(py,j+ph_L)
               if sqrt((i-ss)^2+(s-j)^2)<=ph_L
                  deltapsiL=2;%trucated discrete Laplace with a scale parameter deltaL
                  psi1((s-1)*px+ss,(j-1)*px+i)= 1/(2*deltapsiL)*exp(-sqrt((i-ss)^2+(s-j)^2)/deltapsiL);
               end 
           end
       end  
       for ss=max(1,i-ph_A):min(px,i+ph_A)
           for s=max(1,j-ph_A):min(py,j+ph_A)
               if sqrt((i-ss)^2+(s-j)^2)<=ph_A
                  deltapsiA=2;%trucated discrete Laplace with a scale parameter deltaA
                  psi2((s-1)*px+ss,(j-1)*px+i)= 1/(2*deltapsiA)*exp(-sqrt((i-ss)^2+(s-j)^2)/deltapsiA);
               end    
           end
       end 
   end
end

psi_L=(psi1 ./sum(psi1));
psi_N=psi_L;
psi_A=(psi2 ./sum(psi2));


TE=zeros(patch);
TE1=zeros(patch);
TL=zeros(patch);
TL1=zeros(patch);
TN=zeros(patch);
TN1=zeros(patch);
TA=zeros(patch);
d_A=zeros(patch);





aaa=zeros(k,tlimit);

aaa(:,1)=nzero;
for t=1:tlimit
    jj=mod(t,12);
    if jj==0
        jj=12;
    end
    
    for i=1:patch
    TE(i,i)=S_E*(1-gammaE);  
    TE1(i,i)=S_E*(gammaE);
    for j=1:patch
        TL(i,j)=S_L*m_L(i,j)*(1-gammaL(jj)*sigmaL(jj));
      
        TL1(i,j)=conversion*S_L*psi_L(i,j)*gammaL(jj)*sigmaL(jj);
        TN(i,j)=S_N*m_N(i,j)*(1-gammaN(jj)*sigmaN(jj));
       
        TN1(i,j)=conversion*S_N*psi_N(i,j)*gammaN(jj)*sigmaN(jj);
        TA(i,j)=S_A*m_A(i,j)*(1-gammaA(jj)*sigmaA(jj));
      
        d_A(i,j)=conversion*psi_A(i,j)*gammaA(jj)*sigmaA(jj);
  end
end
  z=zeros(patch);
Th=zeros(patch*4);
TT=zeros(1,patch*4);
Th=[TE z z z;TE1 TL z z; z TL1 TN z;z z TN1 TA];
for i = 1:patch*4
  TT(i) = 1-sum(Th(:,i));
end
Ttilde=[Th; TT];  
    
    
    ff=zeros(k,k);
    for j = 1:patch
           for q = 1:patch  
           ff(j, 3*patch+q) = ((beta(jj)*(1+a*gammaA(jj)*sigmaA(jj)*aaa(3*patch+q,t)))/(1+c*gammaA(j)*sigmaA(jj)*aaa(3*patch+q,t)+a*(gammaA(jj)*sigmaA(jj)*aaa(3*patch+q,t))^2))*d_A(j,q);
           end % q
    end % j
    aaa(:,t+1)=(Th+ff)*aaa(:,t);    
  
end

w=200;
hh=zeros(w,tlimit);
tn=zeros(patch,w,tlimit,'single');
for y=1:w
   fprintf('Current initial number is %d\n',y)
   F=zeros(k,k);
  
n= nzero;
nprime=zeros(k,1);
tn(1,y,1)=sum(nzero');
for itime = 2:tlimit
      jj=mod(itime,12);
    if jj==0
        jj=12;
    end
    
    for i=1:patch
    TE(i,i)=S_E*(1-gammaE);  
    TE1(i,i)=S_E*(gammaE);
    for j=1:patch
        TL(i,j)=S_L*m_L(i,j)*(1-gammaL(jj)*sigmaL(jj));
      
        TL1(i,j)=conversion*S_L*psi_L(i,j)*gammaL(jj)*sigmaL(jj);
        TN(i,j)=S_N*m_N(i,j)*(1-gammaN(jj)*sigmaN(jj));
       
        TN1(i,j)=conversion*S_N*psi_N(i,j)*gammaN(jj)*sigmaN(jj);
        TA(i,j)=S_A*m_A(i,j)*(1-gammaA(jj)*sigmaA(jj));
      
        d_A(i,j)=conversion*psi_A(i,j)*gammaA(jj)*sigmaA(jj);
  end
end
z=zeros(patch);
%Th=zeros(patch*4);
TT=zeros(1,patch*4);
Th=[TE z z z;TE1 TL z z; z TL1 TN z;z z TN1 TA];
for i = 1:patch*4
  TT(i) = 1-sum(Th(:,i));
end
Ttilde=[Th; TT];
    
    nprime=zeros(k,1);
   
    for i = 1:k
     trans = mnrnd(n(i), Ttilde(:, i))';
       nprime = nprime + trans(1:k, :);
       
    for j = 1:patch
           for q = 1:patch  
           F(j, 3*patch+q) = ((beta(jj)*(1+a*gammaA(jj)*sigmaA(jj)*n(3*patch+q)))/(1+c*gammaA(jj)*sigmaA(jj)*n(3*patch+q)+a*(gammaA(jj)*sigmaA(jj)*n(3*patch+q))^2))*d_A(j,q);
           end %q
          if (n(i)>0 && F(j, i)>0)
                  nprime(j) = nprime(j) +poissrnd(F(j, i)*n(i));
            end %if
     
     end %j
     
    end %i
 
     
    n = nprime;
 
hh(y,itime)=sum(nprime(2*patch+1:4*patch));
end

end




hh(:,1)=sum(nzero(2:4)');
figure(3)
plot(sum(aaa(2*patch+1:4*patch,:)),'-k','LineWidth',2)
hold on
plot(mean(hh),'-.k','LineWidth',2)
xlabel('Time')
ylabel('Density of L+N+A')
legend('Deterministic','Stochastic')

Carrying_capacity=mean(sum(aaa(2:4,end-11:end)))
xticks(0:30:240);
xlim([0 240])
ylim([0 170])
function psi_A15=distance_mat2(delta,ph_A,px,py)
counter1=1;
patch=px*py;
for i = 1:px
    %fprintf('Current simulation(psi) number is %d\n',i)
   for j = 1:py
       for ss=max(1,i-ph_A):min(px,i+ph_A)
           for s=max(1,j-ph_A):min(py,j+ph_A)
                tmp=sqrt((i-ss)^2+(s-j)^2);
               if tmp<=ph_A
               
                 counter1=counter1+1;
               end 
           end
       end  
   end
end

psi_A15=zeros(counter1-1,14,'single');
norm_psiA=zeros(patch,1,'single');

counter2=1;
tmp2=1/(2*delta);

%for jj=1:12
for i = 1:px
    
   for j = 1:py
       
       for ss=max(1,i-ph_A):min(px,i+ph_A)
           for s=max(1,j-ph_A):min(py,j+ph_A)
               tmp=sqrt((i-ss)^2+(s-j)^2);
               if tmp<=ph_A
                  row=(s-1)*px+ss;
                  col=(j-1)*px+i;
                  psi_A15(counter2,1)=row;
                  psi_A15(counter2,2)=col;
                  tmp3=tmp2*exp(-tmp/delta);
                  psi_A15(counter2,3)=tmp3;
                  norm_psiA(col)=norm_psiA(col)+ tmp3;
                  counter2=counter2+1;
               end 
           end
       end  
   end
end

for i=1:length(psi_A15(:,1))
  psi_A15(i,3)=psi_A15(i,3)/norm_psiA(psi_A15(i,1));
end
psi_A15=sortrows(psi_A15);
%end
%psi_A15=sortrows(psi_A15);
end
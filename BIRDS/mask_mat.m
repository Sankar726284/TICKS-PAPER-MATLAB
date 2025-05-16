function mask=mask_mat(m_A15)
mask=[1];
   previous=1;
   kb=1;
  for i=1:length(m_A15(:,1))
      current=m_A15(i,1);
      if current~=previous
           %fprintf('Current thnhtjnvtujr is %d\n',current) 
           mask(kb+1)=i;
           kb=kb+1;
      end
      previous=current;
      
  end
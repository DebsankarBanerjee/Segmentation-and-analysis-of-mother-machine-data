function [x1,x2] = cutoff_pos(y,p) 

%%%%%%%%%%%%%%% calculate index(i) where the array 'y(i)' value >= threshold (p)  
%%%%%%%%%%%%%%% from begining x1, from end x2
   
   ny=numel(y);

   for i=1:ny
      if(y(i)>=p)
      x1=i;
      break
      end
   end

   for i=ny:-1:1
      if(y(i)>=p)
      x2=i;
      break
      end
   end


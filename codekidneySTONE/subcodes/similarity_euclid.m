function [R, dmax]= similarity_euclid(data,varargin)


nrow = size(data,1);
R=zeros(nrow,nrow);
data = data';
dmax=0;

if nargin == 1

for i=1:nrow-1
 x=data(:,i);
 for j=i+1:nrow
   y=x-data(:,j);
   d=y'*y;
   d=sqrt(d);
   R(i,j) = d;
   R(j,i) = d; 
   if d>dmax 
       dmax=d; 
   end
 end
end

else

for i=1:nrow-1
 x=data(:,i);
 for j=i+1:nrow
   y=x-data(:,j);
   d=y'*y;
   R(i,j) = d;
   R(j,i) = d; 
   if d>dmax 
       dmax=d; 
   end
 end
end

end
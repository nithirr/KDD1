function R = similarity_pearson(data)

[nrow,ncol] = size(data);
x = mean(data);
data = data-repmat(x,nrow,1);
R = ones(ncol,ncol);

for i = 1:ncol-1
 x = data(:,i);
 X = sqrt(x'*x);
 for j = i+1:ncol
   y = data(:,j);
   xy = x'*y;
   Y = sqrt(y'*y);
   sim = xy/(X*Y);
   R(i,j) = sim;
   R(j,i) = sim; 
 end
end
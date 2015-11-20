function D = distancematrix(loc1,loc2)

if(nargin<2); loc2 = loc1; end

n2 = size(loc2,1);
n1 = size(loc1,1);
a = reshape(loc2,1,n2,2);
b = reshape(loc1,n1,1,2);
D = sqrt(sum((a(ones(n1,1),:,:) - b(:,ones(n2,1),:)).^2,3));


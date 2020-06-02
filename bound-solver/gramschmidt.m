function U = gramschmidt(V)
% GRAMSCHMIDT Creates an orthonormal set of vectors from the input vectors
% using the Gram-Schmidt procedure
%
%   U = GRAMSCHMIDT(V) Creates orthnormal vectors U (as columns) from
%   linearly independent vectors V (as columns)

n = size(V,1);
k = size(V,2);
U = zeros(n,k);
U(:,1) = V(:,1)/sqrt(V(:,1)'*V(:,1));
for nn = 2:k
    U(:,nn) = V(:,nn);
    for mm = 1:nn-1
        U(:,nn) = U(:,nn)-(U(:,mm)'*U(:,nn))/(U(:,mm)'*U(:,mm))*U(:,mm);
    end
    U(:,nn) = U(:,nn)/sqrt(U(:,nn)'*U(:,nn));
end



end
function U = gramschmidt(V)

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
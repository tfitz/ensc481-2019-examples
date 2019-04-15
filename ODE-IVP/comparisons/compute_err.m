function err = compute_err( y0, y1, y1tilde, Atol, Rtol)

n= length(y0);
sc = zeros(n,1);
for i = 1:n
    sc(i) = Atol + max( abs(y0(i)), abs(y1(i)) )*Rtol
end

%% 
err = 1/sqrt(n)*norm( (y1-y1tilde)./sc, 2);
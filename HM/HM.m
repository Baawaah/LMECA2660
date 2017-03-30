% LMECA2660
% Testing zone
fct_cos = @(x) x;%cos(x);
N = 100;
X = linspace(0,1,N);
F_hat = fft(fct_cos(X));
F_hat = F_hat /N;
Y= zeros(1,length(X));
for j = 1 : length(X)
    Y(j) = real(F_hat(1));
    for k = 1 : length(F_hat)/2-1
        Y(j) = 2*( real(F_hat(k))*cos(k*X(j)) - imag(F_hat(k))*sin(k*X(j)) ) + Y(j);  
    end
    Y(j) = 2*( real(F_hat(length(F_hat)/2))*cos( length(F_hat)/2*X(j) ) ) + Y(j);
end
plot(X,Y,X,cos(X));




Q = 100;
H = 5 ;

U1 = @(y) 6*Q/(H*H) * ( y - y.^2 / H );
U2 = @(y) 0.75 * Q /(H/2)* (1 - (y/(H/2)).^2 );

X1 = linspace(0,H,100);
Y1 = U1(X1);

X2 = linspace(-H/2,H/2,100);
Y2 = U2(X2);

plot(X1,Y1,X2,Y2);
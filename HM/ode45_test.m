%Main

A = load('file');
[T,Y] = ode45(@ode45fct,[0 10],[0 0]);
[T Y(:,1)];
A1 = linspace(0,10,100);
A2 = exp(A1)-1;
plot(T,Y(:,1),A(:,1),A(:,2),A1,A2);

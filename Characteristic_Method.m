M = 2.6
y = 1.14
n = 39
nu = 0.5*(sqrt((y+1)/(y-1))*atan(sqrt((y-1)/(y+1)*(M^2-1))) - atan(sqrt(M^2-1)))
theta_max = nu
dtheta = theta_max/n
theta = [dtheta:dtheta:theta_max]
nu = theta
Kn = theta + nu
Kp = theta - nu

Theta = 0.5*bsxfun(@plus,Kn,-Kn')
Nu = 0.5*bsxfun(@minus,Kn,-Kn')
Theta = [theta;Theta]
Nu = [nu;Nu]
Mach = zeros(n+1,n)

for i = 1:n
    for j = 1:(n+1)
        Mach(j,i)=fzero(@(x)(sqrt((y+1)/(y-1))*atan(sqrt((y-1)/(y+1)*(x^2-1)))-atan(sqrt(x^2-1))) - Nu(j,i),2);
    end
end
Mach
Mu = asin(1./Mach)
alpha_p = Theta+Mu
alpha_n = Theta-Mu
dydx_p=tan(alpha_p)
dydx_n=tan(alpha_n)
x = zeros(n+1,n)
y = zeros(n+1,n)

for i=1:n
    for j=i:n
        if i == 1
            if j == 1
                x(i,j)= -2/(dydx_n(i,j)-dydx_p(i,j));
                y(i,j)= 1 + dydx_n(i,j)*x(i,j);
            else
                x(i,j)= (-2)/(dydx_n(i,j)-dydx_p(i,j));
                y(i,j)= 1 + dydx_n(i,j)*x(i,j);
            end
        else
            if i == j
                x(i,j) = x(i-1,j) - y(i-1,j)/dydx_n(i,j);
            else
                x(i,j) = (y(i,j-1) - y(i-1,j) - x(i,j-1)*dydx_p(i,j) + x(i-1,j)*dydx_n(i,j))/(dydx_n(i,j)-dydx_p(i,j));
                y(i,j) = y(i-1,j) + dydx_n(i,j)*(x(i,j)-x(i-1,j));
            end
        end
    end
end

x
y
theta_av = 0.5*(Theta(1:(n),n) + Theta(2:(n+1),n))
a_w = tan(theta_av)
a_c = dydx_p(1:n,n)
X_c = x(1:n,n)
Y_c = y(1:n,n)
X_w = zeros(n+1,1)
Y_w = [1;zeros(n,1)]
for i=2:(n+1)
    X_w(i,1) = (Y_w(i-1,1)-Y_c(i-1,1)-a_w(i-1,1)*X_w(i-1,1)+a_c(i-1,1)*X_c(i-1,1))/(a_c(i-1,1)-a_w(i-1,1));
    Y_w(i,1) = Y_w(i-1,1) + a_w(i-1,1)*(X_w(i,1)-X_w(i-1,1));
end 
X_w
X_c
X_w=X_w(1:(n),1)
Y_w=Y_w(1:(n),1)
plot(X_w,Y_w)
hold on
for i = 3:n
plot([x(i,i:n),X_w(i)],[y(i,i:n),Y_w(i)])
hold on
end
hold off
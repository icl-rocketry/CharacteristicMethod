M = 2.6
y = 1.14
n = 10
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

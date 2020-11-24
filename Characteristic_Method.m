M = 2.6
y = 1.14
n = 10
nu = sqrt((y+1)/(y-1))*atan((y-1)/(y+1)*(M^2-1)) - atan(M^2-1)
theta_max = nu/2
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
        Mach(j,i)=fzero(@(x)sqrt((y+1)/(y-1))*atan(sqrt((y-1)/(y+1))*(x^2-1))-atan(x^2-1) - Nu(j,i),1);
    end
end

Mu = asin(1./Mach)

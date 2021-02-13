%% Housekeeping and initialisation

clear %delete workspace
clc %clear command window
close all %close all figures and other popups

order = 0; %conditions for first iteration
pressurechange = 0; %conditions for first iteration

%READ:
%THIS CODE IS AN IMPLEMENTATION OF THE INVERSE MARCHING METHOD DESCRIBED BY
%ZUCROW AND MODELLED ON THE SHEFFIELD 2020 PAPER WHICH ASSUMES A CENTRELINE
%PRESSURE DISTRIBUTION. 

%THE INVERSE MARCHING METHOD STARTS FROM THE CENTRE LINE AND WORKS BOTTOM
%TO TOP HENCE USING POINTS I,J-1 AND I+1,J TO FIND I,J. THIS METHOD REQUIRES A BUNCH OF ITERATION WHICH IS CAUSING SOME
%NUMERICAL STABILITY PROBLEMS.

%ALSO CURRENTLY THE CODE DOES NOT HAVE GOOD BOUNDARY CONDITIONS AT THE THROAT-
%IT IS NECESSARY TO MODEL A TRANSONIC REGION USING SAUER'S METHOD- SEE THE
%JONES TETT PAPER PAGE 3

%THERE IS AN ALTERNATIVE IN THE DIRECT METHOD WHICH WILL NEED SOME MORE
%BOUNDARY CONDITIONS BUT IT WILL NEED LESS ITERATION SO I WILL LOOK AT THAT
%NEXT

%% Inputs

Lnet = 0.1; %initial guess of centreline "net" length
q = -3; %cubic pressure distribution parameter- must be negative
n = 39; %number of characteristic lines to plot
maxerror = 0.01;

rt = 0.0088; %throat radius

g = 1.135; %ratio of specific heats (gamma)
R = 329;  %gas constant
T0 = 3170; %stagnation temperature
P0 = 30e5; %stagnation pressure
Pe = 1e5; %ambient (and exit) pressure
Cp = R*g/(g-1); %specific heat

Me = sqrt(((P0/Pe)^((g-1)/g)-1)*2/(g-1)); %exit mach number
AR = ((g+1)/2)^((-(g+1))/(2*(g-1)))*((1+(g-1)/2*Me^2)^((g+1)/(2*(g-1))))/Me; %1D area ratio


%% Iteration loop

x = NaN(n,n); %predefine matrices
r = NaN(n,n);
mu = NaN(n,n);
theta = NaN(n,n);
V = NaN(n,n);
T = NaN(n,n);
M = NaN(n,n);
Pn = NaN(n,n);
Qn = NaN(n,n);
N = NaN(n,n);
C = NaN(n,n);

approve = 'y';
while approve == 'y' %iterations
    
    xtemp = x; %save values from last iteration
    rtemp = r;
    Vtemp = V;
    Ttemp = T;
    mutemp = mu;
    thetatemp=theta;
    Mtemp = M;
    
    close all %close all figures
    
    if pressurechange ~= 0
        Lnet = Lnet-Lnet*pressurechange/25; %just shift it into the right direction...
    end
    
    for i=1:n %due to lack of boundary conditions for initial expansion a cubic pressure distribution is assumed along the centreline
        xcentre(i) = Lnet/n*i; %centreline x points
        %         Pcentre(i) = exp(((q*Lnet+2*(log(P0)-log(Pe)))/Lnet^3)*xcentre(i)^3 ...
        %             -((2*q*Lnet + 3*(log(P0)-log(Pe)))/(Lnet^2))*xcentre(i)^2+q*xcentre(i)+log(P0)); %cubic distribution of pressure
        Pcentre(i) = exp(Lnet^(-2)*(xcentre(i)-Lnet)^2*(log(P0)-log(Pe))+log(Pe)); %quadratic distribution of pressure
        if Pcentre(i) >= P0/(1+(g-1)/2)^(g/(g-1)) %if the pressure distribution creates subsonic region
            Pcentre(i) = P0/(1+(g-1)/2*(1+i*0.01)^2)^(g/(g-1)); %make it arbitrarily supersonic instead
        end
        Tcentre(i) = (P0/Pcentre(i))^((1-g)/g)*T0; %centreline properties
        Mcentre(i) = sqrt(((P0/Pcentre(i))^((g-1)/g)-1)*2/(g-1));
        Vcentre(i) = Mcentre(i)*sqrt(g*R*Tcentre(i));
        mucentre(i) = asin(1/Mcentre(i));
        thetacentre(i)=0;
    end
    
    for i= 1:n %for any pair of incident and reflected lines
        x(i,i) = xcentre(i); %pair meets at centreline
        r(i,i) = 0;
        theta(i,i) = thetacentre(i);
        mu(i,i) = mucentre(i);
        M(i,i) = Mcentre(i);
        T(i,i) = Tcentre(i);
        V(i,i) = Vcentre(i);
    end
    
    for k = 1:n %gradually increasing steps away from centreline
        for i = 1:n %for every reflected characteristic
            for j = i+k %for incident characteristics intersecting with reflected ones, stepping away from centreline
                if j >= n+1 %if reached the right end of matrix (crossed all incident characteristics for this reflected characteristic)
                    j = n; %keep j at end
                    break %move to next i
                end
                
                if abs(i-j) == 1 %if just above centreline, then equations are a bit more implicit according to Shapiro 1953 and symbolic toolbox is needed
                    if order == 1 %if iterating with higher order difference
                        
                        Pn(i,j) = theta(i,j)/r(i,j); %arbitrary functions for calculating V
                        Qn(i,j) = theta(i,j)/r(i,j); %Shapiro 1953 shows this is the limit of P and Q as r1, theta1 -> 0 i.e centreline points
                        
                        x(i,j) = (r(i,j-1)-r(i+1,j)+x(i+1,j)*(tan(theta(i+1,j)-mu(i+1,j))+tan(theta(i,j)-mu(i,j)))/2-x(i,j-1)*(tan(theta(i,j-1)+mu(i,j-1))+tan(theta(i,j)+mu(i,j)))/2)... %numerator
                            /((tan(theta(i+1,j)-mu(i+1,j))+tan(theta(i,j)-mu(i,j)))/2-(tan(theta(i,j-1)+mu(i,j-1))+tan(theta(i,j)+mu(i,j)))/2); %denominator
                        
                        r(i,j) = r(i,j-1)+(x(i,j)-x(i,j-1))*(tan(theta(i,j-1)+mu(i,j-1))+tan(theta(i,j)+mu(i,j)))/2; %radial coordinate
                        
                        if j == n %if last incident wave, there is no j+1 point so j value used
                            N(i,j) = 2/(tan(mu(i,j-1))+tan(mu(i,j)))*(2*V(i,j-1)/(V(i,j-1)+V(i,j))+(Pn(i,j)+Pn(i,j))/2*(x(i,j)-x(i,j-1))); %arbitrary functions for calculating V
                            C(i,j) = 2/(tan(mu(i+1,j))+tan(mu(i,j)))*(2*V(i+1,j)/(V(i+1,j)+V(i,j))+(Qn(i,j)+Qn(i,j))/2*(x(i,j)-x(i+1,j))); %because j is at the end the value of j+1 equals value of j.
                            theta(i,j) = theta(i,j-1)+2/(tan(mu(i,j-1))+tan(mu(i,j)))*(2*(V(i,j)-V(i,j-1))/(V(i,j)+V(i,j-1))-(Pn(i,j)+Pn(i,j))/2*(x(i,j)-x(i,j-1))); %flow angle at this point
                        else %all other waves-> have a j+1 point
                            N(i,j) = 2/(tan(mu(i,j-1))+tan(mu(i,j)))*(2*V(i,j-1)/(V(i,j-1)+V(i,j))+(Pn(i,j)+Pn(i,j+1))/2*(x(i,j)-x(i,j-1))); %arbitrary functions for calculating V
                            C(i,j) = 2/(tan(mu(i+1,j))+tan(mu(i,j)))*(2*V(i+1,j)/(V(i+1,j)+V(i,j))+(Qn(i,j)+Qn(i,j+1))/2*(x(i,j)-x(i+1,j))); %arbitrary functions for calculating V
                            theta(i,j) = theta(i,j-1)+2/(tan(mu(i,j-1))+tan(mu(i,j)))*(2*(V(i,j)-V(i,j-1))/(V(i,j)+V(i,j-1))-(Pn(i,j)+Pn(i,j+1))/2*(x(i,j)-x(i,j-1))); %flow angle at this point
                        end
                        
                        V(i,j) = (N(i,j)+C(i,j)-theta(i,j-1)+theta(i+1,j))/(4*((V(i,j-1)+V(i,j))^-1*(tan(mu(i,j-1))+tan(mu(i,j)))^-1+(V(i+1,j)+V(i,j))^-1*(tan(mu(i+1,j))+tan(mu(i,j)))^-1)); %velocity at this point
                        
                    else %if not iterating difference-> first order difference
                        x(i,j) = (r(i,j-1)-r(i+1,j)+x(i+1,j)*tan(theta(i+1,j)-mu(i+1,j))-x(i,j-1)*tan(theta(i,j-1)+mu(i,j-1)))... %x coord, numerator
                            /(tan(theta(i+1,j)-mu(i+1,j))-tan(theta(i,j-1)+mu(i,j-1))); %denominator
                        
                        r(i,j) = r(i,j-1)+(x(i,j)-x(i,j-1))*tan(theta(i,j-1)+mu(i,j-1)); %r coord
                        
                        V(i,j) = (V(i,j-1)*V(i+1,j)*(x(i,j-1)*cot(mu(i,j-1))^2 - x(i,j)*cot(mu(i,j-1))^2 - r(i,j)*theta(i,j-1) + r(i,j)*theta(i+1,j) + theta(i,j-1)*x(i,j-1) - theta(i+1,j)*x(i,j-1) - theta(i,j-1)*x(i,j) + theta(i+1,j)*x(i,j) + r(i,j)*cot(mu(i,j-1)) + r(i,j)*cot(mu(i+1,j)) - x(i,j-1)*cot(mu(i,j-1)) - x(i,j-1)*cot(mu(i+1,j)) + x(i,j)*cot(mu(i,j-1)) + x(i,j)*cot(mu(i+1,j)) - theta(i,j-1)*x(i,j-1)*cot(mu(i,j-1)) - theta(i,j-1)*x(i+1,j)*cot(mu(i+1,j)) + theta(i,j-1)*x(i,j)*cot(mu(i,j-1)) + theta(i,j-1)*x(i,j)*cot(mu(i+1,j)) + x(i+1,j)*cot(mu(i,j-1))*cot(mu(i+1,j)) - x(i,j)*cot(mu(i,j-1))*cot(mu(i+1,j))))... %numerator, this was derived symbolically using P = th(i,j)/r(i,j)
                            /(V(i+1,j)*x(i,j-1)*cot(mu(i,j-1))^2 - V(i+1,j)*x(i,j)*cot(mu(i,j-1))^2 + V(i,j-1)*r(i,j)*cot(mu(i+1,j)) + V(i+1,j)*r(i,j)*cot(mu(i,j-1)) - V(i,j-1)*x(i,j-1)*cot(mu(i+1,j)) - V(i+1,j)*x(i,j-1)*cot(mu(i,j-1)) + V(i,j-1)*x(i,j)*cot(mu(i+1,j)) + V(i+1,j)*x(i,j)*cot(mu(i,j-1)) + V(i+1,j)*x(i+1,j)*cot(mu(i,j-1))*cot(mu(i+1,j)) - V(i+1,j)*x(i,j)*cot(mu(i,j-1))*cot(mu(i+1,j))); %denominator, this was derived symbolically using P = th(i,j)/r(i,j)
                        
                        theta(i,j) = (r(i,j)*(V(i+1,j)*cot(mu(i,j-1))*cot(mu(i+1,j)) - V(i,j-1)*cot(mu(i,j-1))*cot(mu(i+1,j)) + V(i,j-1)*theta(i,j-1)*cot(mu(i+1,j)) + V(i+1,j)*theta(i+1,j)*cot(mu(i,j-1))))... %numerator, this was derived symbolically using P = th(i,j)/r(i,j)
                            /(V(i+1,j)*x(i,j-1)*cot(mu(i,j-1))^2 - V(i+1,j)*x(i,j)*cot(mu(i,j-1))^2 + V(i,j-1)*r(i,j)*cot(mu(i+1,j)) + V(i+1,j)*r(i,j)*cot(mu(i,j-1)) - V(i,j-1)*x(i,j-1)*cot(mu(i+1,j)) - V(i+1,j)*x(i,j-1)*cot(mu(i,j-1)) + V(i,j-1)*x(i,j)*cot(mu(i+1,j)) + V(i+1,j)*x(i,j)*cot(mu(i,j-1)) + V(i+1,j)*x(i+1,j)*cot(mu(i,j-1))*cot(mu(i+1,j)) - V(i+1,j)*x(i,j)*cot(mu(i,j-1))*cot(mu(i+1,j))); %denominator, this was derived symbolically using P = th(i,j)/r(i,j)
                        
                    end
                    T(i,j) = T0 - V(i,j)^2/2/Cp; %temperature
                    M(i,j) = V(i,j)/sqrt(g*R*T(i,j)); %mach number
                    mu(i,j) = asin(1/M(i,j)); %mach angle
                    
                elseif i ~= j %if not on centreline and not just above it, then equations are more explicit and simpler
                    
                    Pn(i,j) = tan(mu(i,j-1))*sin(mu(i,j-1))*sin(theta(i,j-1))/(r(i,j-1)*cos(theta(i,j-1)+mu(i,j-1))); %arbitrary functions for calculating V
                    Qn(i,j) = tan(mu(i+1,j))*sin(mu(i+1,j))*sin(theta(i+1,j))/(r(i+1,j)*cos(theta(i+1,j)-mu(i+1,j))); %arbitrary functions for calculating V
                    
                    if order == 1 %if iterating, use averaged slopes
                        x(i,j) = (r(i,j-1)-r(i+1,j)+x(i+1,j)*(tan(theta(i+1,j)-mu(i+1,j))+tan(theta(i,j)-mu(i,j)))/2-x(i,j-1)*(tan(theta(i,j-1)+mu(i,j-1))+tan(theta(i,j)+mu(i,j)))/2)... %numerator
                            /((tan(theta(i+1,j)-mu(i+1,j))+tan(theta(i,j)-mu(i,j)))/2-(tan(theta(i,j-1)+mu(i,j-1))+tan(theta(i,j)+mu(i,j)))/2); %denominator
                        
                        r(i,j) = r(i,j-1)+(x(i,j)-x(i,j-1))*(tan(theta(i,j-1)+mu(i,j-1))+tan(theta(i,j)+mu(i,j)))/2; %radial component
                        if j == n %if at the end then...
                            N(i,j) = 2/(tan(mu(i,j-1))+tan(mu(i,j)))*(2*V(i,j-1)/(V(i,j-1)+V(i,j))+(Pn(i,j)+Pn(i,j))/2*(x(i,j)-x(i,j-1))); %as before
                            C(i,j) = 2/(tan(mu(i+1,j))+tan(mu(i,j)))*(2*V(i+1,j)/(V(i+1,j)+V(i,j))+(Qn(i,j)+Qn(i,j))/2*(x(i,j)-x(i+1,j)));
                            theta(i,j) = theta(i,j-1)+2/(tan(mu(i,j-1))+tan(mu(i,j)))*(2*(V(i,j)-V(i,j-1))/(V(i,j)+V(i,j-1))-(Pn(i,j)+Pn(i,j))/2*(x(i,j)-x(i,j-1)));
                        else %everywhere else, a j+1 point exists
                            N(i,j) = 2/(tan(mu(i,j-1))+tan(mu(i,j)))*(2*V(i,j-1)/(V(i,j-1)+V(i,j))+(Pn(i,j)+Pn(i,j+1))/2*(x(i,j)-x(i,j-1)));
                            C(i,j) = 2/(tan(mu(i+1,j))+tan(mu(i,j)))*(2*V(i+1,j)/(V(i+1,j)+V(i,j))+(Qn(i,j)+Qn(i,j+1))/2*(x(i,j)-x(i+1,j)));
                            theta(i,j) = theta(i,j-1)+2/(tan(mu(i,j-1))+tan(mu(i,j)))*(2*(V(i,j)-V(i,j-1))/(V(i,j)+V(i,j-1))-(Pn(i,j)+Pn(i,j+1))/2*(x(i,j)-x(i,j-1)));
                        end
                        V(i,j) = (N(i,j)+C(i,j)-theta(i,j-1)+theta(i+1,j))/(4*((V(i,j-1)+V(i,j))^-1*(tan(mu(i,j-1))+tan(mu(i,j)))^-1+(V(i+1,j)+V(i,j))^-1*(tan(mu(i+1,j))+tan(mu(i,j)))^-1));
                    else %not iterating, first order slope
                        x(i,j) = (r(i,j-1)-r(i+1,j)+x(i+1,j)*tan(theta(i+1,j)-mu(i+1,j))-x(i,j-1)*tan(theta(i,j-1)+mu(i,j-1)))/(tan(theta(i+1,j)-mu(i+1,j))-tan(theta(i,j-1)+mu(i,j-1)));
                        r(i,j) = r(i,j-1)+(x(i,j)-x(i,j-1))*tan(theta(i,j-1)+mu(i,j-1));
                        V(i,j) = 1/(cot(mu(i,j-1))/V(i,j-1)+cot(mu(i+1,j))/V(i+1,j))*((cot(mu(i,j-1))*(1+Pn(i,j)*(x(i,j)-x(i,j-1))))+cot(mu(i+1,j))*(1+Qn(i,j)*(x(i,j)-x(i+1,j)))+theta(i+1,j)-theta(i,j-1));
                        theta(i,j) = theta(i,j-1) + cot(mu(i,j-1))*((V(i,j)-V(i,j-1))/V(i,j-1)) - Pn(i,j)*(x(i,j)-x(i,j-1));
                    end
                    
                    T(i,j) = T0 - V(i,j)^2/2/Cp;
                    M(i,j) = V(i,j)/sqrt(g*R*T(i,j));
                    mu(i,j) = asin(1/M(i,j));
                    
                    
                end
                
            end
        end
    end
    
    
    %% Plot characteristics so far
    
    figure() %create figure
    hold on %hold all plots on same figure
    
    for i = 1:n %for all characteristics
        plot(x(i,i:n),r(i,i:n),'r') %plot the reflected plus characteristics
        
        xplot = x; %see below
        xplot(xplot==0)=nan; %data set without ugly 0 values
        
        plot([0,xplot(1,i)],[rt,r(1,i)],'b') %plot the throat-originating minus characteristics before meeting reflections
        plot(xplot(1:n,i),r(1:n,i),'b'); %plot the throat-originating minus characteristics after meeting reflections
        
    end
    
    %% Calculate wall coordinates
    theta_av = zeros(n,1);
    for i = 1:n
        if i <= n-1
            theta_av(i,1) = 0.5*(theta(i,n) + theta(i+1,n)); %flow angle at walls;
        else
            theta_av(i,1) = (theta(i,n));
        end
    end
    %note that n = 39 and we have n+1 lines
    %Theta(:,n) are flow angles after crossing the last throat originating
    %characteristic; these and other values are then constant until they hit
    %the wall. So these are wall angles and values.
    %This line finds the wall flow angle by averaging between two consecutive wall points
    a_w = tan(theta_av); % wall slopes in dy/dx
    a_c = tan(theta(1:n,n)+mu(1:n,n)); % dy/dx slopes of the reflected plus characteristics at wall
    X_c = x(1:n,n); %x coordinate of reflected plus characteristic after crossing last throat-originating minus characteristic
    Y_c = r(1:n,n); %y coordinate of reflected plus characteristic after crossing last throat-originating minus characteristic
    X_w = zeros(n+1,1); %empty matrix of wall x coordinates
    Y_w = [rt;zeros(n,1)]; %empty matrix of wall y coordinates
    
    %% Plot wall
    
    for i=2:n+1 %for all wall points -> start at 2 to n+1 because of numerical scheme propagating from previous point
        %     X_w(i,1) = (Y_w(i-1,1)-Y_c(i-1,1)-a_w(i-1,1)*X_w(i-1,1)+a_c(i-1,1)*X_c(i-1,1))/(a_c(i-1,1)-a_w(i-1,1)); %x coordinate of wall
        %by combining previous wall coordinate and characteristic line
        %meeting the current wall point
        X_w(i,1) = (X_w(i-1,1)*a_w(i-1,1)-X_c(i-1,1)*a_c(i-1,1)+Y_c(i-1,1)-Y_w(i-1,1))/(a_w(i-1,1)-a_c(i-1,1)); %x coordinate of wall
        
        %     Y_w(i,1) = Y_w(i-1,1) + a_w(i-1,1)*(X_w(i,1)-X_w(i-1,1)); %wall y coordinates
        Y_w(i,1) = Y_w(i-1,1) + a_w(i-1,1)*(X_w(i,1)-X_w(i-1,1)); %wall y coordinates
    end
    plot(X_w,Y_w,'k') %plot nozzle wall
    
    for i = 2:n+1
        plot([X_c(i-1),X_w(i)],[Y_c(i-1),Y_w(i)],'r') %plot reflected plus characteristics to wall
    end
    
    if abs(Y_w(end)^2/rt^2 - AR) > maxerror
        disp(['Exit area error= ',num2str((abs(Y_w(end)^2/rt^2 - AR))*100),'%.'])
        approve = input(['Unacceptable discrepancy in exit area, iterate? y/n '],'s');
        approve = string(approve);
        if approve ~= 'y'
            break
        end
        decision = input(['Higher order Difference or new Pressure distribution? d/p '],'s');
        decision = string(decision);
        if decision == 'd'
            order = 1; %iterate slopes
            pressurechange = 0; %do not change centreline distribution
        end
        if decision == 'p'
            pressurechange = Y_w(end)^2/rt^2 - AR; %the error in area ratio used to change centreline distribution
            order = 0; %do not iterate slopes
        end
        
        
    else disp(['discrepancy between calculated and 1-D area ratio less than ',num2str(maxerror*100),'%, analysis complete'])
        approve = 'n'; %do not iterate
    end
end

%% Postprocessing

plot([0,X_w(end)],[0,0],'k','linestyle','- -') %plot centreline
pointsize = 50; %scatter point size
scatter(x(:),r(:),pointsize,M(:)) %create scatter heatmap for mach number
cb = colorbar; %create colorbar for mach scatter plot
title(cb,'Mach number')
title(['Axisymmetric nozzle. Exit area error of ',num2str(abs(Y_w(end)^2/rt^2 - AR)*100),'% '])

errorx = xtemp - x; %change in x from last iteration
errorxreal = real(errorx); %real change
errorximag = imag(errorx); %imaginary change (non physical, the sim is screwed)

figure()
plot(1:n,errorxreal)
title('real residuals in x for each reflected line (last iteration)')
leg = legend;
title(leg,'index of reflected line')
figure()
plot(1:n,errorximag)
title('imaginary components (if not 0 then the current sim is screwed)')


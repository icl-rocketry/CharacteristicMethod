%% Housekeeping

clear %delete workspace
clc %clear command window
close all %close all figures and other popups

%% Definitions and initialisation

M = 2.6; %Exit Mach number (obtained from engine pressure ratio and ambient pressure, assuming perfet expansion)
y = 1.14; %ratio of specific heats (gamma)
n = 39; %number of characteristic lines to plot
nu = 0.5*(sqrt((y+1)/(y-1))*atan(sqrt((y-1)/(y+1)*(M^2-1))) - atan(sqrt(M^2-1))); %P-M function of half of exit mach number
theta_max = nu; %set max angle to half the P-M function of exit mach number (this is the maximum nozzle turning angle and the angle of the initial corner)

dtheta = theta_max/n; %steps in angles
theta = [dtheta:dtheta:theta_max]; %array of angles
nu = theta; %set nu to be array of angles
Kn = theta + nu; %constant for right running wave (top left to bottom right)
Kp = theta - nu; %constant for left running wave (bottom left to top right)

Theta = 0.5*bsxfun(@plus,Kn,-Kn'); %theta obtained back from Kn at every intersection point
Nu = 0.5*bsxfun(@minus,Kn,-Kn'); %nu obtained back from Kn at every intersection point
    %bsxfun essentially intermeshes the arrays into matrices
    
Theta = [theta;Theta]; %matrix of theta (flow angles) at each point
Nu = [nu;Nu]; %matrix of nu (PM function angles) at each point
Mach = zeros(n+1,n); %empty matrix placeholder for mach numbers

for i = 1:n %for all discrete angle steps
    for j = 1:(n+1) %for all discrete angle steps plus 1
        Mach(j,i)=fzero(@(x)(sqrt((y+1)/(y-1))*atan(sqrt((y-1)/(y+1)*(x^2-1)))-atan(sqrt(x^2-1))) - Nu(j,i),2); % use fzero function 
%         to find M solution of PM function for given nu, so this is Mach number at every intersection point. 
    end
end

Mu = asin(1./Mach); %matrix of mach angles (physically the angle in which information propagates (mach cone for example) and the angle of expansion fans)
alpha_p = Theta+Mu; %constants for plus characteristic (bottom left to top right)
alpha_n = Theta-Mu; %constants for minus characteristic (top left to bottom right)
dydx_p=tan(alpha_p); %slopes of plus characteristic line
dydx_n=tan(alpha_n); %slopes of minus characteristic line 
x = zeros(n+1,n); %empty matrix of x positions -> these will be points of intersecting characteristics that plot the lines
y = zeros(n+1,n); %empty matrix of y positons -> these will be points of intersecting characteristics that plot the lines

%% Calculate points on reflected plus characteristics

%i denotes the index of throat-originated minus characteristics
%j denotes the index of the reflected plus characteristics

for i=1:n %for all reflected plus characteristics
    for j=i:n %for all throat-originated minus characteristics that this reflected line crosses
        if i == 1 %if first reflected plus characteristic
            if j == 1 %if first throat-originated minus characteristic (has some numerical funkiness)
                x(i,j)= -2/(dydx_n(i,j)-dydx_p(i,j)); %dx after moving dy = -1 (i.e throat to centreline) 
                y(i,j)= 1 + dydx_n(i,j)*x(i,j); %going from 1 down to centreline with slope; should always equal 0.
            else
                x(i,j)= (-2)/(dydx_n(i,j)-dydx_p(i,j)); %dx after moving dy = -1 (i.e throat to centreline)
                y(i,j)= 1 + dydx_n(i,j)*x(i,j); %going from 1 down to centreline with slope; should always equal 0.
            end
        else %for all other reflected plus characteristics after first
            if i == j %point of reflection i.e on centreline
                x(i,j) = x(i-1,j) - y(i-1,j)/dydx_n(i,j); %calculate from previous point and current line slope knowing that y = 0 at centreline.
            else %points of intersections with throat-originating minus characteristics
                x(i,j) = (y(i,j-1) - y(i-1,j) - x(i,j-1)*dydx_p(i,j) + x(i-1,j)*dydx_n(i,j))/(dydx_n(i,j)-dydx_p(i,j)); %see below 
                   %This the exact same equation used in method of characteristics
                   %to find the location of points from two preceding points from
                   %which two characteristic lines originate and which intersect at
                   %the point i,j. Simply connect the dots.
                y(i,j) = y(i-1,j) + dydx_n(i,j)*(x(i,j)-x(i-1,j)); %knowing x(i,j) you only need to use one previous point to find y.
            end
        end
    end
end


%% Plot characteristics so far

figure() %create figure
hold on %hold all plots on same figure

for i = 2:n %for all characteristics
plot(x(i,i:n),y(i,i:n),'r') %plot the reflected plus characteristics

xplot = x; %see below
xplot(xplot==0)=nan; %data set without ugly 0 values

plot([0,xplot(1,i)],[1,y(1,i)],'b') %plot the throat-originating minus characteristics before meeting reflections
plot(xplot(1:n,i),y(1:n,i),'b'); %plot the throat-originating minus characteristics after meeting reflections

end

%% Calculate wall coordinates

theta_av = 0.5*(Theta(1:(n),n) + Theta(2:(n+1),n)); %flow angle at walls;
    %note that n = 39 and we have n+1 lines
    %Theta(:,n) are flow angles after crossing the last throat originating
    %characteristic; these and other values are then constant until they hit
    %the wall. So these are wall angles and values.
    %This line finds the wall flow angle by averaging between two consecutive wall points
a_w = tan(theta_av); % wall slopes in dy/dx
a_c = dydx_p(1:n,n); % dy/dx slopes of the reflected plus characteristics at wall
X_c = x(1:n,n); %x coordinate of reflected plus characteristic after crossing last throat-originating minus characteristic
Y_c = y(1:n,n); %y coordinate of reflected plus characteristic after crossing last throat-originating minus characteristic
X_w = zeros(n+1,1); %empty matrix of wall x coordinates
Y_w = [1;zeros(n,1)]; %empty matrix of wall y coordinates

%% Plot wall

for i=2:(n+1) %for all wall points -> start at 2 to n+1 because of numerical scheme propagating from previous point
    X_w(i,1) = (Y_w(i-1,1)-Y_c(i-1,1)-a_w(i-1,1)*X_w(i-1,1)+a_c(i-1,1)*X_c(i-1,1))/(a_c(i-1,1)-a_w(i-1,1)); %x coordinate of wall
        %by combining previous wall coordinate and characteristic line
        %meeting the current wall point
    Y_w(i,1) = Y_w(i-1,1) + a_w(i-1,1)*(X_w(i,1)-X_w(i-1,1)); %wall y coordinates
end 
plot(X_w,Y_w,'k') %plot nozzle wall
for i = 2:n+1
plot([X_c(i-1),X_w(i)],[Y_c(i-1),Y_w(i)],'r') %plot reflected plus characteristics to wall
end

%% Postprocessing

plot([0,X_w(end)],[0,0],'k','linestyle','- -') %plot centreline

pointsize=15;
scatter(x(:),y(:),pointsize,Mach(:))
colorbar
% throatradius = 3.4;
% Xgeometry = X_w*throatradius;
% Ygeometry = Y_w*throatradius;

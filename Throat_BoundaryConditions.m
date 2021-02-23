%Scrip to work ou tthe boundary conditions at the throat of the nozzle
%using the Keiran's paper

%%
clear
clc
%%
n=20;
R_t = 0.088;
rho_t = R_t*2;
gamma = 1.14;
T = 3000;
R = 329;
a_s = sqrt((2*gamma*R*T)/(gamma + 1));
P = 3000000;

r = linspace(0.0, 0.088, 30);

a = sqrt(2/((gamma+1)*(rho_t*R_t)));
E = -((R_t)/(8))*((2*(gamma+1)/((rho_t/R_t))))^0.5;

U = ((gamma+1)*(a^2).*r)/2;
x_s = -((gamma+1)*a.*r.^2)/(4);
x_i = -((gamma+1)*a.*r.^2)/(8);

x_s_E = x_s - E;
x_i_E = x_i - E;



u = a.*x_i + ((gamma+1)/4)*a^2 * r.^2;
v = (gamma+1)/2 * a^2 .*x_i .*r + (((gamma+1)^2)/16)*a^3 .*r.^3;

u_bar = a_s*(1+u);
v_bar = a_s * v;

M = u_bar/a_s;


t = T./(1+((gamma-1)/(2)).*M.^2);
p = P./(1+((gamma-1)/(2)).*M.^2).^(gamma/(gamma-1));
rho = p./(R.*t);
%%
a_0 = gamma*329*3000;

figure
%plot(u, v)
xlabel("u")
ylabel("v")


figure

non_dim_height = r./R_t;
non_dim_length_i = x_i_E./R_t;
non_dim_length_s = x_s./R_t;
plot(non_dim_length_i, non_dim_height, "LineWidth", 1)
hold on
% plot(non_dim_length_s, non_dim_height)
xlabel("Non dimensional length x/R_t")
ylabel("Non dimensional length r/R_t")
xlim([0 0.3])
ylim([0 1])
% figure
% plot(u, r)


mu(:,1) = asin(1./M');
drdx_p = tan(mu);
drdx_m = tan(-mu);

% y = mx+c

for i=1:length(r)
    x1 = non_dim_length_i(i);
    y1 = non_dim_height(i);
    
    x2 = x1 + 0.01;
    y2 = y1 + 0.01*drdx_p(i);
    y3 = y1 + 0.01*drdx_m(i);
    
    p1 = polyfit([x1, x2], [y1, y2], 1);
    c = p1(2);
    f = polyval(p1, [x1, x2]);
    plot([x1, x2], f, "-r", "LineWidth", 1.4);
    
    p2 = polyfit([x1, x2], [y1, y3], 1);
    c_1 = p2(2);
    f1 = polyval(p2, [x1, x2]);
    plot([x1, x2], f1, "-r", "LineWidth", 1.4)
end
theta = zeros(n,n);
%%

for i = 2:n %for every reflected characteristic
    for j = i:n %for incident characteristics intersecting with reflected ones, stepping away from centreline
        if i==j
            Pn(i,j) = tan(mu(i,j-1))*sin(mu(i,j-1))*sin(theta(i,j-1))/(r(i,j-1)*cos(theta(i,j-1)+mu(i,j-1))); %arbitrary functions for calculating V
            %             Qn(i,j) = tan(mu(i-1,j))*sin(mu(i-1,j))*sin(theta(i-1,j))/(r(i-1,j)*cos(theta(i-1,j)-mu(i-1,j))); %arbitrary functions for calculating V
            
            
            x(i,j) = (r(i,j-1)-x(i,j-1)*tan(theta(i,j-1)+mu(i,j-1)))/(-tan(theta(i,j-1)+mu(i,j-1)));
            r(i,j) = r(i,j-1)+(x(i,j)-x(i,j-1))*tan(theta(i,j-1)+mu(i,j-1));
            
            V(i,j) = V(i-1,j)/(1+(theta(i-1,j)-theta(i,j))*tan(mu(i,j))+((tan(mu(i,j))*sin(mu(i,j))*sin(theta(i,j))*(x(i-1,j)-x(i,j)))/(cos(theta(i,j)-mu(i,j))*r(i,j))));
            
            
%             V(i,j) = 1/(cot(mu(i,j-1))/V(i,j-1)+cot(mu(i-1,j))/V(i-1,j))*((cot(mu(i,j-1))*(1+Pn(i,j)*(x(i,j)-x(i,j-1))))+cot(mu(i-1,j))*(1+Qn(i,j)*(x(i,j)-x(i-1,j)))+theta(i-1,j)-theta(i,j-1));
            theta(i,j) = theta(i,j-1) + cot(mu(i,j-1))*((V(i,j)-V(i,j-1))/V(i,j-1)) - Pn(i,j)*(x(i,j)-x(i,j-1));
        end
        
        
        Pn(i,j) = tan(mu(i,j-1))*sin(mu(i,j-1))*sin(theta(i,j-1))/(r(i,j-1)*cos(theta(i,j-1)+mu(i,j-1))); %arbitrary functions for calculating V
        Qn(i,j) = tan(mu(i-1,j))*sin(mu(i-1,j))*sin(theta(i-1,j))/(r(i-1,j)*cos(theta(i-1,j)-mu(i-1,j))); %arbitrary functions for calculating V
        
        
        x(i,j) = (r(i,j-1)-r(i-1,j)+x(i-1,j)*tan(theta(i-1,j)-mu(i-1,j))-x(i,j-1)*tan(theta(i,j-1)+mu(i,j-1)))/(tan(theta(i-1,j)-mu(i-1,j))-tan(theta(i,j-1)+mu(i,j-1)));
        r(i,j) = r(i,j-1)+(x(i,j)-x(i,j-1))*tan(theta(i,j-1)+mu(i,j-1));
        
        V(i,j) = 1/(cot(mu(i,j-1))/V(i,j-1)+cot(mu(i-1,j))/V(i-1,j))*((cot(mu(i,j-1))*(1+Pn(i,j)*(x(i,j)-x(i,j-1))))+cot(mu(i-1,j))*(1+Qn(i,j)*(x(i,j)-x(i-1,j)))+theta(i-1,j)-theta(i,j-1));
        theta(i,j) = theta(i,j-1) + cot(mu(i,j-1))*((V(i,j)-V(i,j-1))/V(i,j-1)) - Pn(i,j)*(x(i,j)-x(i,j-1));
        
        T(i,j) = T0 - V(i,j)^2/2/Cp;
        M(i,j) = V(i,j)/sqrt(g*R*T(i,j));
        mu(i,j) = asin(1/M(i,j));
    end
end


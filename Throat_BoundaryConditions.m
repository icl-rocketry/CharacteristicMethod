%Scrip to work ou tthe boundary conditions at the throat of the nozzle
%using the Keiran's paper

%%
clear
clc
close all

%%
n=20;
R_t = 0.088;
rho_t = R_t*2;
gamma = 1.14;
T = 3000;
R = 329;
a_s = sqrt((2*gamma*R*T)/(gamma + 1));
P = 30e5;
T0 = 3000;
r = linspace(0.0, 0.088, n);

a = sqrt(2/((gamma+1)*(rho_t*R_t)));
E = -((R_t)/(8))*((2*(gamma+1)/((rho_t/R_t))))^0.5;

U = ((gamma+1)*(a^2).*r)/2; %
x_s = -((gamma+1)*a.*r.^2)/(4); %supersonic line
x_i = -((gamma+1)*a.*r.^2)/(8); %initial value line

x_s_E = x_s - E;
x_i_E = x_i - E;

u = a.*x_i + ((gamma+1)/4)*a^2 * r.^2; %nondimensionalised
v = (gamma+1)/2 * a^2 .*x_i .*r + (((gamma+1)^2)/16)*a^3 .*r.^3;

u_bar = a_s*(1+u); %u_bar is the physical velocity in x
v_bar = a_s * v; %in y

M = u_bar/a_s;


t = T./(1+((gamma-1)/(2)).*M.^2);
p = P./(1+((gamma-1)/(2)).*M.^2).^(gamma/(gamma-1));
rho = p./(R.*t);
Cp = R*gamma/(gamma-1); %specific heat
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


% mu(:,1) = asin(1./M');
% drdx_p = tan(mu);
% drdx_m = tan(-mu);

% y = mx+c;
%
% for i=1:length(r)
%     x1 = non_dim_length_i(i);
%     y1 = non_dim_height(i);
%
%     x2 = x1 + 0.01;
%     y2 = y1 + 0.01*drdx_p(i);
%     y3 = y1 + 0.01*drdx_m(i);
%
%     p1 = polyfit([x1, x2], [y1, y2], 1);
%     c = p1(2);
%     f = polyval(p1, [x1, x2]);
%     plot([x1, x2], f, "-r", "LineWidth", 1.4);
%
%     p2 = polyfit([x1, x2], [y1, y3], 1);
%     c_1 = p2(2);
%     f1 = polyval(p2, [x1, x2]);
%     plot([x1, x2], f1, "-r", "LineWidth", 1.4)
% end
x_i  = x_i_E;
r_i  = r;
mu_i = fliplr(asin(1./M));
theta_i = zeros(1,n);
% V_i = fliplr(sqrt((1+u).^2+v.^2))*a_o;
T_i = fliplr(t);
M_i = fliplr(M);
V_i = fliplr(M.*sqrt(gamma*R.*t));

x = NaN(2*n,n); %predefine matrices
r = NaN(2*n,n);
mu = NaN(2*n,n);
theta = NaN(2*n,n);
V = NaN(2*n,n);
T = NaN(2*n,n);
M = NaN(2*n,n);


for i = 1:n %populate initial value line
    x(i,i) =  x_i(i);
    r(i,i) = r_i(i);
    mu(i,i) = mu_i(i);
    theta(i,i) = theta_i(i);
    V(i,i) = V_i(i);
    T(i,i) = T_i(i);
    M(i,i) = M_i(i);
end

% %%
% i = 1;
% j = 2; %tricky point because i,j-1 has r2 = 0, limit needed



%% %% Initial Value line block
for k = 1:n %for every reflected characteristic
    for i = 2:n %for incident characteristics intersecting with reflected ones, stepping away from centreline
        if i + k <= n %if j does not exceed n
            j = i+k;
            if i~=j %if not on the initial line
                
                r2 = r(i+1,j);
                r1 = r(i,j-1);
                x2 = x(i+1,j);
                x1 = x(i,j-1);
                mu2 = mu(i+1,j);
                mu1 = mu(i,j-1);
                th2 = theta(i+1,j);
                th1 = theta(i,j-1);
                V2 = V(i+1,j);
                V1 = V(i,j-1);
                
                %x3 = (r2 - r1 + x1*tan(mu1 + th1) + x2*tan(mu2 - th2))/(tan(mu2 - th2) + tan(mu1 + th1))
                x(i,j) = (r2 - r1 + x1*tan(mu1 + th1) + x2*tan(mu2 - th2))/(tan(mu2 - th2) + tan(mu1 + th1));
                x3 = x(i,j);
                %r3 = r1 - x1*tan(mu1 + th1) + x3*tan(mu1 + th1)
                r(i,j) = r1 - x1*tan(mu1 + th1) + x3*tan(mu1 + th1);
                r3 = r(i,j);
                
                if i == 1 && j == 2
                    
                    r1 = r3;
                    r2 = r3;
                    
                end
                
                
                V(i,j) = (V1*r1*cos(mu1 + th1) - V1*r1*th1*cos(mu1 + th1)*tan(mu1) - V1*x1*sin(mu1)*tan(mu1)*sin(th1) + V1*x3*sin(mu1)*tan(mu1)*sin(th1))/(r1*cos(mu1 + th1)) - (V1*(r1*r3*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu1) + r1*r3*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu2) - r1*r3*th1*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu1)*tan(mu2) - r1*r3*th2*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu1)*tan(mu2) + r1*x2*cos(mu1 + th1)*sin(mu2)*tan(mu1)*tan(mu2)*sin(th2) - r1*x3*cos(mu1 + th1)*sin(mu2)*tan(mu1)*tan(mu2)*sin(th2) + r1*r3*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu1)*tan(mu2)*((V1*r1*r3*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu1) + V1*r1*r3*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu2) - V1*r3*x1*cos(mu2 - th2)*sin(mu1)*tan(mu1)*tan(mu2)*sin(th1) + V1*r3*x3*cos(mu2 - th2)*sin(mu1)*tan(mu1)*tan(mu2)*sin(th1) - V1*r1*r3*th1*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu1)*tan(mu2) - V1*r1*r3*th2*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu1)*tan(mu2) + V1*r1*x2*cos(mu1 + th1)*sin(mu2)*tan(mu1)*tan(mu2)*sin(th2) - V1*r1*x3*cos(mu1 + th1)*sin(mu2)*tan(mu1)*tan(mu2)*sin(th2))^2/(V1^2*r1^2*r3^2*cos(mu2 - th2)^2*cos(mu1 + th1)^2*tan(mu1)^2*tan(mu2)^2) + (4*(V2*r1*r3*cos(mu2 - th2)*cos(mu1 + th1) - V1*r1*r3*cos(mu2 - th2)*cos(mu1 + th1) + V1*r3*x1*cos(mu2 - th2)*sin(mu1)*tan(mu1)*sin(th1) - V1*r3*x3*cos(mu2 - th2)*sin(mu1)*tan(mu1)*sin(th1) + V1*r1*r3*th1*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu1) + V1*r1*r3*th2*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu2) - V1*r1*x2*cos(mu1 + th1)*sin(mu2)*tan(mu2)*sin(th2) + V1*r1*x3*cos(mu1 + th1)*sin(mu2)*tan(mu2)*sin(th2) + V1*x3^2*sin(mu1)*sin(mu2)*tan(mu1)*tan(mu2)*sin(th1)*sin(th2) + V1*r1*th1*x2*cos(mu1 + th1)*sin(mu2)*tan(mu1)*tan(mu2)*sin(th2) - V1*r1*th1*x3*cos(mu1 + th1)*sin(mu2)*tan(mu1)*tan(mu2)*sin(th2) - V1*r3*th2*x1*cos(mu2 - th2)*sin(mu1)*tan(mu1)*tan(mu2)*sin(th1) + V1*r3*th2*x3*cos(mu2 - th2)*sin(mu1)*tan(mu1)*tan(mu2)*sin(th1) - V1*r1*r3*th1*th2*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu1)*tan(mu2) + V1*x1*x2*sin(mu1)*sin(mu2)*tan(mu1)*tan(mu2)*sin(th1)*sin(th2) - V1*x1*x3*sin(mu1)*sin(mu2)*tan(mu1)*tan(mu2)*sin(th1)*sin(th2) - V1*x2*x3*sin(mu1)*sin(mu2)*tan(mu1)*tan(mu2)*sin(th1)*sin(th2)))/(V1*r1*r3*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu1)*tan(mu2)))^(1/2) - r3*x1*cos(mu2 - th2)*sin(mu1)*tan(mu1)*tan(mu2)*sin(th1) + r3*x3*cos(mu2 - th2)*sin(mu1)*tan(mu1)*tan(mu2)*sin(th1)))/(2*r1*r3*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu2));                
                V3 = V(i,j);
                theta(i,j) = ((V2 - V3)/V3 + th2*tan(mu2) - (sin(mu2)*tan(mu2)*sin(th2)*(x2 - x3))/(r3*cos(mu2 - th2)))/tan(mu2);                T(i,j) = T0 - V(i,j)^2/2/Cp;
                M(i,j) = V(i,j)/sqrt(gamma*R*T(i,j));
                mu(i,j) = asin(1/M(i,j));
            end
        end
    end
    %values(n+1,j) = values(1,j)
end
%% Centreline block

for i = n+1:2*n
    for k =1:n
        if i+k <2*n
            j = i-n+k;
            if i-n ==j %point on centreline
                %
            else %not on centreline
                r2 = r(i-1,j);
                r1 = r(i,j-1);
                x2 = x(i-1,j);
                x1 = x(i,j-1);
                mu2 = mu(i-1,j);
                mu1 = mu(i,j-1);
                th2 = theta(i-1,j);
                th1 = theta(i,j-1);
                V2 = V(i-1,j);
                V1 = V(i,j-1);
                
                
                %x3 = (r2 - r1 + x1*tan(mu1 + th1) + x2*tan(mu2 - th2))/(tan(mu2 - th2) + tan(mu1 + th1))
                x(i,j) = (r2 - r1 + x1*tan(mu1 + th1) + x2*tan(mu2 - th2))/(tan(mu2 - th2) + tan(mu1 + th1));
                x3 = x(i,j);
                %r3 = r1 - x1*tan(mu1 + th1) + x3*tan(mu1 + th1)
                r(i,j) = r1 - x1*tan(mu1 + th1) + x3*tan(mu1 + th1);
                r3 = r(i,j);
                
                
                % V3 = (V1*r1*cos(mu1 + th1) + V1*r1*th1*cos(mu1 + th1)*tan(mu1) - V1*x1*sin(mu1)*tan(mu1)*sin(th1) + V1*x3*sin(mu1)*tan(mu1)*sin(th1))/(r1*cos(mu1 + th1)) - (V1*(r1*r2*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu1) + r1*r2*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu2) + r1*r2*th1*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu1)*tan(mu2) + r1*r2*th2*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu1)*tan(mu2) - r1*x2*cos(mu1 + th1)*sin(mu2)*tan(mu1)*tan(mu2)*sin(th2) + r1*x3*cos(mu1 + th1)*sin(mu2)*tan(mu1)*tan(mu2)*sin(th2) + r1*r2*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu1)*tan(mu2)*((V1*r1*r2*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu1) + V1*r1*r2*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu2) - V1*r2*x1*cos(mu2 - th2)*sin(mu1)*tan(mu1)*tan(mu2)*sin(th1) + V1*r2*x3*cos(mu2 - th2)*sin(mu1)*tan(mu1)*tan(mu2)*sin(th1) + V1*r1*r2*th1*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu1)*tan(mu2) + V1*r1*r2*th2*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu1)*tan(mu2) - V1*r1*x2*cos(mu1 + th1)*sin(mu2)*tan(mu1)*tan(mu2)*sin(th2) + V1*r1*x3*cos(mu1 + th1)*sin(mu2)*tan(mu1)*tan(mu2)*sin(th2))^2/(V1^2*r1^2*r2^2*cos(mu2 - th2)^2*cos(mu1 + th1)^2*tan(mu1)^2*tan(mu2)^2) - (4*(V1*r1*r2*cos(mu2 - th2)*cos(mu1 + th1) - V2*r1*r2*cos(mu2 - th2)*cos(mu1 + th1) - V1*r2*x1*cos(mu2 - th2)*sin(mu1)*tan(mu1)*sin(th1) + V1*r2*x3*cos(mu2 - th2)*sin(mu1)*tan(mu1)*sin(th1) + V1*r1*r2*th1*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu1) + V1*r1*r2*th2*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu2) - V1*r1*x2*cos(mu1 + th1)*sin(mu2)*tan(mu2)*sin(th2) + V1*r1*x3*cos(mu1 + th1)*sin(mu2)*tan(mu2)*sin(th2) + V1*x3^2*sin(mu1)*sin(mu2)*tan(mu1)*tan(mu2)*sin(th1)*sin(th2) - V1*r1*th1*x2*cos(mu1 + th1)*sin(mu2)*tan(mu1)*tan(mu2)*sin(th2) + V1*r1*th1*x3*cos(mu1 + th1)*sin(mu2)*tan(mu1)*tan(mu2)*sin(th2) - V1*r2*th2*x1*cos(mu2 - th2)*sin(mu1)*tan(mu1)*tan(mu2)*sin(th1) + V1*r2*th2*x3*cos(mu2 - th2)*sin(mu1)*tan(mu1)*tan(mu2)*sin(th1) + V1*r1*r2*th1*th2*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu1)*tan(mu2) + V1*x1*x2*sin(mu1)*sin(mu2)*tan(mu1)*tan(mu2)*sin(th1)*sin(th2) - V1*x1*x3*sin(mu1)*sin(mu2)*tan(mu1)*tan(mu2)*sin(th1)*sin(th2) - V1*x2*x3*sin(mu1)*sin(mu2)*tan(mu1)*tan(mu2)*sin(th1)*sin(th2)))/(V1*r1*r2*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu1)*tan(mu2)))^(1/2) - r2*x1*cos(mu2 - th2)*sin(mu1)*tan(mu1)*tan(mu2)*sin(th1) + r2*x3*cos(mu2 - th2)*sin(mu1)*tan(mu1)*tan(mu2)*sin(th1)))/(2*r1*r2*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu2))
                V(i,j) = (V1*r1*cos(mu1 + th1) + V1*r1*th1*cos(mu1 + th1)*tan(mu1) - V1*x1*sin(mu1)*tan(mu1)*sin(th1) + V1*x3*sin(mu1)*tan(mu1)*sin(th1))/(r1*cos(mu1 + th1)) - (V1*(r1*r2*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu1) + r1*r2*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu2) + r1*r2*th1*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu1)*tan(mu2) + r1*r2*th2*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu1)*tan(mu2) - r1*x2*cos(mu1 + th1)*sin(mu2)*tan(mu1)*tan(mu2)*sin(th2) + r1*x3*cos(mu1 + th1)*sin(mu2)*tan(mu1)*tan(mu2)*sin(th2) + r1*r2*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu1)*tan(mu2)*((V1*r1*r2*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu1) + V1*r1*r2*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu2) - V1*r2*x1*cos(mu2 - th2)*sin(mu1)*tan(mu1)*tan(mu2)*sin(th1) + V1*r2*x3*cos(mu2 - th2)*sin(mu1)*tan(mu1)*tan(mu2)*sin(th1) + V1*r1*r2*th1*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu1)*tan(mu2) + V1*r1*r2*th2*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu1)*tan(mu2) - V1*r1*x2*cos(mu1 + th1)*sin(mu2)*tan(mu1)*tan(mu2)*sin(th2) + V1*r1*x3*cos(mu1 + th1)*sin(mu2)*tan(mu1)*tan(mu2)*sin(th2))^2/(V1^2*r1^2*r2^2*cos(mu2 - th2)^2*cos(mu1 + th1)^2*tan(mu1)^2*tan(mu2)^2) - (4*(V1*r1*r2*cos(mu2 - th2)*cos(mu1 + th1) - V2*r1*r2*cos(mu2 - th2)*cos(mu1 + th1) - V1*r2*x1*cos(mu2 - th2)*sin(mu1)*tan(mu1)*sin(th1) + V1*r2*x3*cos(mu2 - th2)*sin(mu1)*tan(mu1)*sin(th1) + V1*r1*r2*th1*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu1) + V1*r1*r2*th2*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu2) - V1*r1*x2*cos(mu1 + th1)*sin(mu2)*tan(mu2)*sin(th2) + V1*r1*x3*cos(mu1 + th1)*sin(mu2)*tan(mu2)*sin(th2) + V1*x3^2*sin(mu1)*sin(mu2)*tan(mu1)*tan(mu2)*sin(th1)*sin(th2) - V1*r1*th1*x2*cos(mu1 + th1)*sin(mu2)*tan(mu1)*tan(mu2)*sin(th2) + V1*r1*th1*x3*cos(mu1 + th1)*sin(mu2)*tan(mu1)*tan(mu2)*sin(th2) - V1*r2*th2*x1*cos(mu2 - th2)*sin(mu1)*tan(mu1)*tan(mu2)*sin(th1) + V1*r2*th2*x3*cos(mu2 - th2)*sin(mu1)*tan(mu1)*tan(mu2)*sin(th1) + V1*r1*r2*th1*th2*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu1)*tan(mu2) + V1*x1*x2*sin(mu1)*sin(mu2)*tan(mu1)*tan(mu2)*sin(th1)*sin(th2) - V1*x1*x3*sin(mu1)*sin(mu2)*tan(mu1)*tan(mu2)*sin(th1)*sin(th2) - V1*x2*x3*sin(mu1)*sin(mu2)*tan(mu1)*tan(mu2)*sin(th1)*sin(th2)))/(V1*r1*r2*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu1)*tan(mu2)))^(1/2) - r2*x1*cos(mu2 - th2)*sin(mu1)*tan(mu1)*tan(mu2)*sin(th1) + r2*x3*cos(mu2 - th2)*sin(mu1)*tan(mu1)*tan(mu2)*sin(th1)))/(2*r1*r2*cos(mu2 - th2)*cos(mu1 + th1)*tan(mu2));
                V3 = V(i,j);
                % th3 = ((V1 - V3)/V1 + th1*tan(mu1) - (sin(mu1)*tan(mu1)*sin(th1)*(x1 - x3))/(r1*cos(mu1 + th1)))/tan(mu1)
                theta(i,j) = ((V1 - V3)/V1 + th1*tan(mu1) - (sin(mu1)*tan(mu1)*sin(th1)*(x1 - x3))/(r1*cos(mu1 + th1)))/tan(mu1);
                
                T(i,j) = T0 - V(i,j)^2/2/Cp;
                M(i,j) = V(i,j)/sqrt(gamma*R*T(i,j));
                mu(i,j) = asin(1/M(i,j));
            end
            
        end
    end
end


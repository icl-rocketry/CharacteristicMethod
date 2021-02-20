%Scrip to work ou tthe boundary conditions at the throat of the nozzle
%using the Keiran's paper

%%
clear
clc
%%

R_t = 0.088;
rho_t = R_t*2;
gamma = 1.14;

r = linspace(0.0, 0.088, 50);

a = sqrt(2/((gamma+1)*(rho_t*R_t)));
E = -((R_t)/(8))*((2*(gamma+1)/((rho_t/R_t))))^0.5;

U = ((gamma+1)*(a^2).*r)/2;
M = U+1;
x_s = -((gamma+1)*a.*r.^2)/(4)-E;
hold on
x_i = -((gamma+1)*a.*r.^2)/(8) - E;
hold off


u = a.*x + ((gamma+1)/4)*a^2 * r.^2;
v = (gamma+1)/2 * a^2 .*x .*r + (((gamma+1)^2)/16)*a^3 .* r.^3;
%%
a_0 = gamma*329*3000;

figure
plot(u, v)
xlabel("u")
ylabel("v")


figure

non_dim_height = r./R_t;
non_dim_length = x./R_t;
plot(non_dim_length, non_dim_height)
xlabel("Non dimensional length x/R_t")
ylabel("Non dimensional length r/R_t")

% figure
% plot(u, r)
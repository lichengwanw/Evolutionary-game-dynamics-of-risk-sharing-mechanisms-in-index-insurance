clear
clc
close all
color_library = [[233, 196, 107];  
    [230, 111, 081]; 
    [038, 070, 083]; 
    [042, 157, 142]]./255;
%% Parameter assignment
w = 1; 
p = 0.2;
q = 0.2;
alpha = 0.7;
r = 0.001;
gamma = 0.8; 
detalist = [0,0.04,0.1,0.2]; 
c = 0.17; 
N = 50; 
beta = 0.9;
M=10;
UW =@(w) w^(1-gamma)/(1-gamma); 

%% The replication equation selects the gradient
syms x
for rr=1:length(detalist)
    deta = detalist(rr);
    combMatrix = NaN(N+1, N+1);  
    for i = 0:N
        for j = 0:min(i, N)
            combMatrix(i+1, j+1) = nchoosek(i, j);
        end
    end

    E_NO_CII = (1-p)*UW(w)+p*UW((1-alpha)*w);
    f_x = 0;
    for k=0:N-1
        f_x = f_x+nchoosek(N-1,k)* x^k * (1-x)^(N-1-k)*pi_C(k+1,alpha,w,c,deta,UW,q,p,r,combMatrix,M,beta);
    end
    f_x = f_x-E_NO_CII;
    f_x = x*(1-x)*f_x;
    x_vals = 0:0.01:1;
    y_vals = subs(f_x, x, x_vals);
    plot(x_vals, y_vals, 'Color', color_library(rr, :), 'LineWidth', 2);
    hold on;
    legendInfo{rr} = ['$\delta = ', num2str(deta), '$'];  
end
plot(x_vals, zeros(1, length(x_vals)), 'k', 'HandleVisibility', 'off'); 
legend(legendInfo, 'Interpreter', 'latex', 'FontSize', 20,'Location', 'best');
hold off;
ax = gca;
ax.FontSize = 25;  
xticks(linspace(0, 1, 6));
ax.YAxis.Exponent = -2;  
pbaspect([1.2 1 1]);

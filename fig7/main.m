clear
clc
close all

color_library = [[233, 196, 107];  
                 [230, 111, 081];  
                 [038, 070, 083];  
                 [042, 157, 142]]./255;
n = 64;
my_colormap = zeros(n, 3);
my_colormap(:, 1) = linspace(1, 0, n);
my_colormap(1:n/2, 2) = linspace(0, 1, n/2);
my_colormap(n/2+1:end, 2) = linspace(1, 0, n/2);
my_colormap(:, 3) = linspace(0, 1, n);
colormap(flipud(my_colormap));


%% Parameter assignment
w = 1; 
p = 0.2;
q = 0.2;
alpha = 0.7;
r = 0.001;
gamma = 0.8; 
detalist = [0,0.05,0.1,0.2]; 
betas = 10;
c = 0.17; 
Z = 50; 
u = 0.02;
N = 40; 
beta = 0.9;
M=10;
Average_insurance_participation_rate = [];
set(0, 'DefaultAxesFontName', 'Times New Roman');

for ff = 1:length(detalist)
    deta = detalist(ff);
    UW =@(w) w^(1-gamma)/(1-gamma); 
    combMatrix = NaN(Z+1, Z+1); 
    for i = 0:Z
        for j = 0:min(i, N)
            combMatrix(i+1, j+1) = nchoosek(i, j);
        end
    end
    TC_N =@(i,Fc,Fn)  (1-u)*(i/(Z))*((Z-i)/(Z-1))*((1+exp(betas*(Fc-Fn)))^(-1))+u*(i/Z);%%C-N
    TN_C =@(i,Fn,Fc)  (1-u)*(i/(Z))*((Z-i)/(Z-1))*((1+exp(betas*(Fn-Fc)))^(-1))+u*((Z-i)/Z);%%N-C

    p_shift = zeros(Z+1,Z+1); 
    E_NO_CII = (1-p)*UW(w)+p*UW((1-alpha)*w);
    len=size(p_shift,1);
    p_diff = zeros(1,len);
    for i=1:len-1
        p_shift(i,i+1) = TN_C(i-1,E_NO_CII,F_C_CII(i-1,Z,N,alpha,w,c,deta,UW,q,p,r,combMatrix,M,beta));
        p_shift(i+1,i) = TC_N(i,F_C_CII(i,Z,N,alpha,w,c,deta,UW,q,p,r,combMatrix,M,beta),E_NO_CII); 
    end
    for i=1:len
        p_diff(i) = TN_C(i-1,E_NO_CII,F_C_CII(i-1,Z,N,alpha,w,c,deta,UW,q,p,r,combMatrix,M,beta))-TC_N(i-1,F_C_CII(i-1,Z,N,alpha,w,c,deta,UW,q,p,r,combMatrix,M,beta),E_NO_CII);
    end
    aaaa = sum(p_shift');
    for i=1:Z+1
        p_shift(i,i) = 1-aaaa(i);
    end
    sss = ones(1,Z+1)*(p_shift-diag(ones(Z+1,1))+ones(Z+1,Z+1))^(-1);
    stage = 0:Z;
    Average_insurance_participation_rate = [Average_insurance_participation_rate;sum(sss.*stage)/Z];

    figure(1)
    plshow = plot(stage, sss, '.', 'MarkerSize', 20, 'Color', color_library(mod(ff-1, size(color_library,1))+1, :), 'DisplayName', ['$\delta = $ ', num2str(deta)]);
    hold on
    for i=1:length(stage)
        line([stage(i) stage(i)], [0 sss(i)], 'LineStyle', '-', 'Color', color_library(mod(ff-1, size(color_library,1))+1, :),'LineWidth', 1,'HandleVisibility', 'off');
    end
    legend('show', 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Times New Roman')
    pbaspect([1.2 1 1]);
    xlim([0,50])
    ylim([0,0.2])
    xticks(linspace(0, 50, 6));
    ax = gca;
    ax.FontSize = 25;

    figure(2)
    plshow2 = plot(stage,p_diff,'.','MarkerSize', 20, 'DisplayName', ['$delta = $', num2str(deta)], 'Color', color_library(mod(ff-1, size(color_library,1))+1, :));
    hold on
    xlim([0,50])
    pbaspect([1.2 1 1]);
    xticks(linspace(0, 50, 6));
    
end

plot(stage,zeros(1,Z+1),'k', 'HandleVisibility','off')  %
ax = gca;
ax.FontSize = 25;
ylim([-0.02,0.08])
yticks(linspace(-0.02,0.08, 5));
ax.YAxis.Exponent = -2;  

figure
b = bar(1:length(detalist),Average_insurance_participation_rate,'FaceColor', 'flat');
b.FaceColor = 'flat';
b.CData = color_library;

formatted_labels = arrayfun(@(x) sprintf('$%.2f$', x), detalist, 'UniformOutput', false);
set(gca, 'TickLabelInterpreter', 'latex'); 
xticklabels(formatted_labels); 
yticks(linspace(0,1,6))
xtips1 = b.XEndPoints;
ytips1 = b.YEndPoints;
labels1 = string(round(b.YData, 2));
text(xtips1, ytips1, labels1, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 20)
pbaspect([1.2 1 1]);
ax = gca;
ax.FontSize = 25;
xlim([0, 5])
ylim([0,1.05])

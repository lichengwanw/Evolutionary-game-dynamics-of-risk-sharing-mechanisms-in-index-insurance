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
rlist = [0.001,0.05,0.1 0.15]; 
gamma = 0.8; 
detalist = linspace(0, 0.12, 22); 
betas = 10; 
c = 0.17; 
Z = 50; 
u = 0.02;
N = 40; 
beta = 0.9; 
M = 10;

%% Markov process
figure;
hold on;
legends = {}; 
for rr = 1:length(rlist)
    r = rlist(rr);
    for ff = 1:length(detalist)
        deta = detalist(ff);
        UW = @(w) w^(1-gamma) / (1-gamma); 
        combMatrix = NaN(Z+1, Z+1); 
        for i = 0:Z
            for j = 0:min(i, N)
                combMatrix(i+1, j+1) = nchoosek(i, j);
            end
        end

        TC_N = @(i, Fc, Fn) (1-u) * (i/(Z)) * ((Z-i)/(Z-1)) * ((1+exp(betas*(Fc-Fn)))^(-1)) + u * (i/Z); %% C-N
        TN_C = @(i, Fn, Fc) (1-u) * (i/(Z)) * ((Z-i)/(Z-1)) * ((1+exp(betas*(Fn-Fc)))^(-1)) + u * ((Z-i)/Z); %% N-C

        p_shift = zeros(Z+1, Z+1); 
        E_NO_CII = (1-p) * UW(w) + p * UW((1-alpha) * w);
        len = size(p_shift, 1);
        p_diff = zeros(1, len); 
        for i = 1:len-1
            p_shift(i, i+1) = TN_C(i-1, E_NO_CII, F_C_CII(i-1, Z, N, alpha, w, c, deta, UW, q, p, r, combMatrix, M, beta));
            p_shift(i+1, i) = TC_N(i, F_C_CII(i, Z, N, alpha, w, c, deta, UW, q, p, r, combMatrix, M, beta), E_NO_CII);
        end
        for i = 1:len
            p_diff(i) = TN_C(i-1, E_NO_CII, F_C_CII(i-1, Z, N, alpha, w, c, deta, UW, q, p, r, combMatrix, M, beta)) - TC_N(i-1, F_C_CII(i-1, Z, N, alpha, w, c, deta, UW, q, p, r, combMatrix, M, beta), E_NO_CII);
        end
        aaaa = sum(p_shift');
        for i = 1:Z+1
            p_shift(i, i) = 1-aaaa(i);
        end
        sss = ones(1, Z+1) * (p_shift - diag(ones(Z+1, 1)) + ones(Z+1, Z+1))^(-1);
        stage = 0:1:Z;
        Average_insurance_participation_rate = 0;
        for i = 1:length(stage)
            Average_insurance_participation_rate = Average_insurance_participation_rate + stage(i) * sss(i);
        end
        z(ff) = Average_insurance_participation_rate/Z;
    end

    if rr == 1
        plot(detalist, z, '-d', 'LineWidth', 2, 'MarkerSize', 6, 'Color', color_library(rr,:));
    elseif rr == 2
        plot(detalist, z, '-*', 'LineWidth', 2, 'MarkerSize', 6, 'Color', color_library(rr,:));
    elseif rr == 3
        plot(detalist, z, '-h', 'LineWidth', 2, 'MarkerSize', 6, 'Color', color_library(rr,:));
    else 
        plot(detalist, z, '-o', 'LineWidth', 2, 'MarkerSize', 6, 'Color', color_library(rr,:));
    end
    legends{end+1} = ['$r = ', num2str(r), '$'];
end

legend(legends, 'Interpreter', 'latex', 'FontSize', 15,'Location', 'best');


hold off;
ax = gca;
ax.FontSize = 17; 
xticks(linspace(min(detalist), max(detalist), 7));
yticks(linspace(0, 1, 6));
ax.Box = 'on';
xlim([min(detalist), max(detalist)]);
ylim([0,1])
pbaspect([1.2 1 1]);


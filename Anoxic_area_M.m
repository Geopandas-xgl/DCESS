load('OutThilda_M')

global dm d n 

% Get parameter values
ParVal_M    % Activate global parameters

zcent= [dm/2 dm+(d:d:(n-1)*d)-d/2]/1e3;   % Vertical center of boxes
%% 
figure

% Define the range for x and y
x = st/1e3;
y = zcent;

% Create the grid
[X, Y] = meshgrid(x, y);

% Compute the function values
Z = permute(sLL(8,:,:),[2,3,1])*1000;

% Create the contour plot
[c, h] = contour(X, Y, Z);
clabel(c, h, 'FontSize', 6)
colorbar;
xlabel('Time (kyr)');
ylabel('Depth (km)');
title("Figure 1");
set(gca,"YDir","reverse");
%% 
GAw = load('GAw_Eo.txt');
area_layer = nan(55,1);
area_layer(1) = 0;
for i = 1:54
   area_layer(i+1) = (GAw(i)-GAw(i+1))/GAw(1);
end

area_layer_matrix = repmat(area_layer,1,length(x));

figure("Position",[0,0,1000,300])
[h,~]= tight_subplot(1,3,[0,0.08],[0.15,0.1],[0.08,0.01]);
axes(h(1))
for i = 20:20:120
area_layer_matrix_1 = area_layer_matrix;
area_layer_matrix_1(Z>=i) = 0; % 0.2 ml/mol: 8.93 mmol/m3; 0.5 ml/mol: 22.3 mmol/m3
area_anox_1 = sum(area_layer_matrix_1)*100;
plot(st/1e3,area_anox_1, "LineWidth",1.2); hold on
xlabel("Time (kyr)");
ylabel("Area fraction (%)");
end
lgd = legend(h(1),string(20:20:120));
lgd.Title.String = "O_{2}";
title("Figure 2a");
set(gca,"linewidth", 1,"FontSize",12,"FontName", "Times New Roman","TickLength",[0.02,0.025],"Layer","top");

axes(h(2))
for i = 20
area_layer_matrix_1 = area_layer_matrix;
area_layer_matrix_1(Z>=i) = 0; % 0.2 ml/mol: 8.93 mmol/m3; 0.5 ml/mol: 22.3 mmol/m3
area_anox_1 = sum(area_layer_matrix_1)*100;
plot(st/1e3,area_anox_1, "LineWidth",1.2); hold on
xlabel("Time (kyr)");
ylabel("Area fraction (%)");
end
lgd = legend(h(2),string(20));
lgd.Title.String = "O_{2}";
title("Figure 2b");
set(gca,"linewidth", 1,"FontSize",12,"FontName", "Times New Roman","TickLength",[0.02,0.025],"Layer","top");

axes(h(3))
oxy = 60:20:80;
clr = ["#F2C965","#7E2F8E"];
for i = 1: length(oxy)
oxy_i = oxy(i); 
area_layer_matrix_1 = area_layer_matrix;
area_layer_matrix_1(Z>=oxy_i) = 0; % 0.2 ml/mol: 8.93 mmol/m3; 0.5 ml/mol: 22.3 mmol/m3
area_anox_1 = sum(area_layer_matrix_1)*100;
fold_anox_1 = area_anox_1/area_anox_1(1);
plot(st/1e3,fold_anox_1,"Color",clr(i),"LineWidth",1.2); hold on
xlabel("Time (kyr)");
ylabel("expanding folds");
end
legend(h(3),string(60:20:80));
title("Figure 2c");
set(gca,"linewidth", 1,"FontSize",12,"FontName", "Times New Roman","TickLength",[0.02,0.025],"Layer","top");






% oxy_layer_matrix = Z.*area_layer_matrix;
% oxy_sum = sum(oxy_layer_matrix(3:30,:));
% oxy_norm = oxy_sum/oxy_sum(1); % In the modern seawater, the fraction of  seafloor anoxia (<0.2 ml/mol) is 0.2% 
% anox_area = 1./oxy_norm;
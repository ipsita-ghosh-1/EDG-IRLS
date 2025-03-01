% r = load('tilde05-Aug-2024 19:39:32_alg_name__type_GaussianData_ovrsmp_1_4_rmax_5')
% r_1 = load('tilde+105-Aug-2024 20:11:38_alg_name__type_GaussianData_ovrsmp_1_4_rmax_5.mat')
% r_2 = load('tilde+205-Aug-2024 20:39:26_alg_name__type_GaussianData_ovrsmp_1_4_rmax_5.mat')
% 
% r_2r = load('tilde2_05-Aug-2024 21:04:42_alg_name__type_GaussianData_ovrsmp_1_4_rmax_5.mat')

% r_3r = load('02-Aug-2024_3r_22:32:53_alg_name__type_GaussianData_ovrsmp_2_4_rmax_5.mat')
error_r = median(r.rel_error_Dist,4);
error_r1 = median(r_1.rel_error_Dist,4);
error_r2 = median(r_2.rel_error_Dist,4);
error_2r = median(r_2r.rel_error_Dist,4);
% error_3r = abs(r_3r.rel_error_Dist);
oversampling_list=[1,1.5,2,2.5,3,3.5,4];
% oversampling_list = [1.0000,   1.5000    2.0000    2.5000    3.0000    3.5000    4.0000]



figure;
% Plot the first line
semilogy(oversampling_list,error_2r(:,:,1),'-p','LineWidth', 2, 'DisplayName', 'MatrixIRLS');
% % Plot the second line
hold on;
semilogy(oversampling_list,error_2r(:,:,2),'-p','LineWidth', 2, 'DisplayName', 'ScaledSGD');
hold on;
semilogy(oversampling_list,error_2r(:,:,3),'-p','LineWidth', 2, 'DisplayName', 'ALM');
hold on; 
semilogy(oversampling_list,error_2r(:,:,4),'-p','LineWidth', 2, 'DisplayName', 'RieEDG');
% hold on;
% semilogy(oversampling_list,error_3r(:,:,4),'-p','LineWidth', 2, 'DisplayName', '3r');
% legend('show');
% Adjust y-axis limits to squeeze the plot
ylim([min(error_2rs(:))/2, max(error_2r(:))*2]);

% Customize legend
lgd = legend('show');
set(lgd, 'FontSize', 18, 'Location', 'best');
% Add labels and title with larger and bolder font
xlabel('Oversampling Rate', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Median Relative Error', 'FontSize', 14, 'FontWeight', 'bold');

% Make the numbers on x and y axes larger and bolder
set(gca, 'FontSize', 22, 'FontWeight', 'bold');
% Add labels and title
% xlabel('oversampling rate');
% ylabel('Median Relative Error')

% 
% figure;
% % Plot the first line
% semilogy(oversampling_list,error_r(:,:,1),'-p','LineWidth', 2, 'DisplayName', 'r');
% % % Plot the second line
% hold on;
% semilogy(oversampling_list,error_r1(:,:,1),'-p','LineWidth', 2, 'DisplayName', 'r+1');
% hold on;
% semilogy(oversampling_list,error_r2(:,:,1),'-p','LineWidth', 2, 'DisplayName', 'r+2');
% hold on; 
% semilogy(oversampling_list,error_2r(:,:,1),'-p','LineWidth', 2, 'DisplayName', '2r');
% legend('show');
% 
% xlabel('oversampling rate');
% ylabel('Relative Error')

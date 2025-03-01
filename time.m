
m= length(outs_cell{1,1}{1,1}.time)
s = length(outs_cell{1,1}{1,2}.time)
a= length(outs_cell{1,1}{1,3}.time)
r = length(outs_cell{1,1}{1,4}.time)
max_length = max([m, s, a, r]);
array_from_1_to_max = 1:max_length;
mlen = 1:m
figure;
% Plot the first line
plot(mlen,outs_cell{1,1}{1,1}.time,'-p','LineWidth', 2, 'DisplayName', 'MatrixIRLS');
% % Plot the second line
hold on;
plot(array_from_1_to_max,outs_cell{1,1}{1,2}.time,'-p','LineWidth', 2, 'DisplayName', 'ScaledSGD');
% % Plot the second line
hold on;
plot(array_from_1_to_max,outs_cell{1,1}{1,3}.time,'-p','LineWidth', 2, 'DisplayName', 'ALM');
% % Plot the second line
hold on;
plot(array_from_1_to_max,outs_cell{1,1}{1,4}.time,'-p','LineWidth', 2, 'DisplayName', 'ReiEDG');
legend('show');
% ylim([0,200])
xlabel('Iterations');
ylabel('Runtime')
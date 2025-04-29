%% 2.1
syms a b;
s = -243*(-1/2*a-1/2*b)-32*(2*a+3*b+1/12)-(-5/2*a-15/2*b-2/3)+243*(1/2*a-1/2*b)+32*(-2*a+3*b-1/12)+(5/2*a-15/2*b+2/3)
s = 729*(-1/2*a-1/2*b)+64*(2*a+3*b+1/12)+(-5/2*a-15/2*b-2/3)+729*(1/2*a-1/2*b)+64*(-2*a+3*b-1/12)+(5/2*a-15/2*b+2/3)
s = -2187*(-1/2*a-1/2*b)-128*(2*a+3*b+1/12)-(-5/2*a-15/2*b-2/3)+2187*(1/2*a-1/2*b)+128*(-2*a+3*b-1/12)+(5/2*a-15/2*b+2/3)

%% 2.2
syms a b k;
kk = (-1/2*a-1/2*b)*(cos(-3*k)+i*sin(-3*k))+(2*a+3*b+1/12)*(cos(-2*k)+i*sin(-2*k))+(-5/2*a-15/2*b-2/3)*(cos(-1*k)+i*sin(-1*k))...
+10*b+(5/2*a-15/2*b+2/3)*(cos(k)+i*sin(k))+(-2*a+3*b-1/12)*(cos(2*k)+i*sin(2*k))+(1/2*a-1/2*b)*(cos(3*k)+i*sin(3*k));
kkk = kk/i;
m = simplify(kkk);
pretty(m)

k = 0:0.01:pi;
a = [-0.2,0,1/30,0.1,0.2];
for i = 1:length(a)
    Re = (4/3+5*a(i))*sin(k) - (1/6+4*a(i))*sin(2*k)+a(i)*sin(3*k);
    plot(k,Re,'LineWidth',1.5);
    hold on;
end
legend('\alpha=-0.2','\alpha=0','\alpha=1/30','\alpha=0.1','\alpha=0.2','Location','best');
xlabel('k');
ylabel('Re(k'')');
title('不同\alpha取值下的色散曲线');

b = [-1,0,0.1,1,2];
for i = 1:length(b)
    Im = 15*b(i)*cos(k)-6*b(i)*cos(2*k)+b(i)*cos(3*k)-10*b(i);
    plot(k,Im,'LineWidth',1.5);
    hold on;
end
legend('\beta=-1','\beta=0','\beta=0.1','\beta=1','\beta=2','Location','best');
xlabel('k');
ylabel('Im(k'')');
title('不同\beta取值下的耗散曲线');

%% 2.3
nu = [4 6 8 10];
a = -0.1:0.00001:0.1;

hold on;
for i = 1:length(nu)
    E = zeros(size(a)); % 预分配结果数组
    % 定义被积函数（确保标量输入输出）
    Ee = @(k, a_val) exp(nu(i)*(pi-k)) .* ((4/3+5*a_val)*sin(k) - (1/6+4*a_val)*sin(2*k) + a_val*sin(3*k) - k).^2;
    
    % 对每个a值计算积分
    for j = 1:length(a)
        E(j) = integral(@(k) Ee(k, a(j)), 0, pi, 'AbsTol', 1e-5) / exp(nu(i)*pi);
    end
    plot(a, log10(E), 'LineWidth', 1.5);
    
    % 找到极小值点
    [min_E(i), idx] = min(E);
    min_a(i) = a(idx);
    
    % 绘制极小值点的虚线
    xline(min_a(i), '--', 'LineWidth', 0.5, 'LabelHorizontalAlignment', 'center');
end
hold off;

xlabel('\alpha');
% 使用 TeX 格式
ytick_labels = {'10^{-10}', '10^{-9}','10^{-8}', '10^{-7}', '10^{-6}', '10^{-5}', '10^{-4}', '10^{-3}'};
set(gca, 'YTickLabel', ytick_labels, 'TickLabelInterpreter', 'tex');
ylabel('E');
legend('\nu=4', '','\nu=6', '','\nu=8', '','\nu=10','','Location','best');
title('目标函数随\alpha变化曲线');
grid on;

fprintf('极小值点:\n');
for i = 1:length(nu)
    fprintf('nu = %d: a = %.5f, E = %.4e\n', nu(i), min_a(i), min_E(i));
end
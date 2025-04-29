%% 3 初始条件
%clear;
Mx=256;%空间网格数
Nt=15000;%时间网格数
N=Mx*Nt; %总网格数
delta_x = 1/Mx;
delta_t = 10/Nt;
CFL = 1*delta_t/delta_x;
c = 1;%波速

U=zeros(1,Mx+1);
%先赋初值和边界值
for x = 1:Mx+1
    U(x)=Initial3(20,0+(x-1)*delta_x);
end
U(Mx+1)=U(1);%边界周期条件

%% DRP格式
a = [-0.02651995, 0.18941314, -0.79926643, 0, 0.79926643, -0.18941314, 0.02651995];
b = [2.30255809, -2.49100760, 1.57434093, -0.38589142];

% 历史右端项存储（四层显式方法需要前三步）
rhs_history = zeros(4, Mx+1); % 存储前三步的右端项

% 前三步使用RK4启动
for t = 1:3
    U1 = U-1/4*c*delta_t*DRPcompute_rhs(U, Mx, delta_x, a);
    U2 = U-1/3*c*delta_t*DRPcompute_rhs(U1, Mx, delta_x, a);
    U3 = U-1/2*c*delta_t*DRPcompute_rhs(U2, Mx, delta_x, a);
    U = U-c*delta_t*DRPcompute_rhs(U3, Mx, delta_x, a);
    rhs_history(t+1, :) = DRPcompute_rhs(U, Mx, delta_x, a); % 存储历史项
end

% 后续使用四层显式方法推进
for t = 4:Nt
    % 计算当前右端项
    rhs_current = DRPcompute_rhs(U, Mx, delta_x, a);
    % 更新历史队列
    rhs_history = circshift(rhs_history, [-1,0]);
    rhs_history(4, :) = rhs_current;

    % 四层显式更新公式
    U = U - delta_t * (b(1)*rhs_history(4,:) + b(2)*rhs_history(3,:) + ...
              b(3)*rhs_history(2,:) + b(4)*rhs_history(1,:));
end
U(Mx+1) = U(1);
DRP = U;


%% DRP-M格式
a = [-0.020843142770, 0.166705904415, -0.770882380518, 0, ...
     0.770882380518, -0.166705904415, 0.020843142770]; % j=-3到3
c_damp = [-0.014281184692, 0.086150669577, -0.235718815308, 0.327698660846, ...
          -0.235718815308, 0.086150669577, -0.014281184692]; % j=-3到3
b = [2.30255809, -2.49100760, 1.57434093, -0.38589142];
nu_a = 0.001*delta_x^2;

% 历史右端项存储（四层显式方法需要前三步）
rhs_history = zeros(4, Mx+1); % 存储前三步的右端项

% 前三步使用RK4启动
for t = 1:3
    U1 = U-1/4*c*delta_t*DRPcompute_rhs(U, Mx, delta_x, a);
    U2 = U-1/3*c*delta_t*DRPcompute_rhs(U1, Mx, delta_x, a);
    U3 = U-1/2*c*delta_t*DRPcompute_rhs(U2, Mx, delta_x, a);
    U = U-c*delta_t*DRPcompute_rhs(U3, Mx, delta_x, a);
    rhs_history(t+1, :) = DRPMcompute_rhs(U, Mx, delta_x, a, c, nu_a, c_damp); % 存储历史项
end

% 后续使用四层显式方法推进
for t = 4:Nt
    % 计算当前右端项
    rhs_current = DRPMcompute_rhs(U, Mx, delta_x, a, c,nu_a, c_damp);
    % 更新历史队列
    rhs_history = circshift(rhs_history, [-1,0]);
    rhs_history(4, :) = rhs_current;

    % 四层显式更新公式
    U = U + delta_t * (b(1)*rhs_history(4,:) + b(2)*rhs_history(3,:) + ...
              b(3)*rhs_history(2,:) + b(4)*rhs_history(1,:));
end
U(Mx+1) = U(1);
DRPM = U;

%% MDCD格式
a = 0.0463783;
b = 0.012;
% 使用RK4
for t = 1:Nt
    U1 = U-1/4*c*delta_t*MDCDcompute_dudx(U, Mx, delta_x, a, b);
    U2 = U-1/3*c*delta_t*MDCDcompute_dudx(U1, Mx, delta_x, a, b);
    U3 = U-1/2*c*delta_t*MDCDcompute_dudx(U2, Mx, delta_x, a, b);
    U = U-c*delta_t*MDCDcompute_dudx(U3, Mx, delta_x, a, b);
end
MDCD = U;

%% SA-DRP格式
epsilon = 1e-08;
% 使用RK4
for t = 1:Nt
    U1 = U-1/4*c*delta_t*SADRPcompute_dfdx(U, Mx, delta_x, epsilon);
    U2 = U-1/3*c*delta_t*SADRPcompute_dfdx(U1, Mx, delta_x, epsilon);
    U3 = U-1/2*c*delta_t*SADRPcompute_dfdx(U2, Mx, delta_x, epsilon);
    U = U-c*delta_t*SADRPcompute_dfdx(U3, Mx, delta_x, epsilon);
end
SADRP = U;



%% 作图
X = 0:delta_x:1;
plot(X,SADRP,'LineWidth',1.0);
xlabel('x','Fontsize',14);
ylabel('u','Fontsize',14);
legend('t=10.0','Fontsize',14);
title('SA-DRP格式数值解','Fontsize',14);
grid on;

%% 综合作图
exact = 0;
subplot(1,2,2);
XX = 0.8:delta_x:1;
X = 0.8:0.001:1;
for i = 1:20
    exact = exact + sin(2*pi*i*(X-10));
end
exact = exact/20;
hold on;
plot(X,exact,'k');
scatter(XX,DRP(206:257),30,'bo');
scatter(XX,DRPM(206:257),30,'rs');
scatter(XX,MDCD(206:257),30,'m^');
scatter(XX,SADRP(206:257),30,'gv');

plot(XX,DRP(206:257),'b','LineWidth',1.0);
plot(XX,DRPM(206:257),'r','LineWidth',1.0);
plot(XX,MDCD(206:257),'m','LineWidth',1.0);
plot(XX,SADRP(206:257),'g','LineWidth',1.0);

xlabel('x','Fontsize',14);
ylabel('u','Fontsize',14);
legend('exact','DRP','DRP-M','MDCD','SA-DRP','Location','best','Fontsize',14);
title('多种半离散格式数值解（部分）','Fontsize',14);
hold off;

subplot(1,2,1);
exact = 0;
XX = 0:delta_x:1;
X = 0:0.001:1;
for i = 1:20
    exact = exact + sin(2*pi*i*(X-10));
end
exact = exact/20;
hold on;
plot(X,exact,'k');
scatter(XX,DRP,10,'bo');
scatter(XX,DRPM,10,'rs');
scatter(XX,MDCD,10,'m^');
scatter(XX,SADRP,10,'gv');

plot(XX,DRP,'b','LineWidth',1.0);
plot(XX,DRPM,'r','LineWidth',1.0);
plot(XX,MDCD,'m','LineWidth',1.0);
plot(XX,SADRP,'g','LineWidth',1.0);

xlabel('x','Fontsize',14);
ylabel('u','Fontsize',14);
legend('exact','DRP','DRP-M','MDCD','SA-DRP','Location','best','Fontsize',14);
title('多种半离散格式数值解','Fontsize',14);
hold off;

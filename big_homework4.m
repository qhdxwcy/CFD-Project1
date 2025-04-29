%% 4 初始条件
%clear;
Mx=256;%空间网格数
Nt=Mx*5;%时间网格数
N=Mx*Nt; %总网格数
delta_x = 1/Mx;
delta_t = 1/Nt;
CFL = 1*delta_t/delta_x
c = 1;%波速

%先单独取一组随机数，保证所有格式初值条件相同，用后注释
% r = zeros(1,64);
% for i = 1:64
%     r(i) = rand;
% end
U=zeros(1,Mx+1);
%先赋初值和边界值
for x = 1:Mx+1
    U(x)=Initial4(0.1,24,r,0+(x-1)*delta_x);
    %U(x)=Initial3(20,0+(x-1)*delta_x);
end
U(Mx+1)=U(1);%边界周期条件
U0 = U;

%% DRP格式
a = [-0.02651995, 0.18941314, -0.79926643, 0, 0.79926643, -0.18941314, 0.02651995];

% 使用RK4
for t = 1:Nt
    U1 = U-1/4*c*delta_t*DRPcompute_rhs(U, Mx, delta_x, a);
    U2 = U-1/3*c*delta_t*DRPcompute_rhs(U1, Mx, delta_x, a);
    U3 = U-1/2*c*delta_t*DRPcompute_rhs(U2, Mx, delta_x, a);
    U = U-c*delta_t*DRPcompute_rhs(U3, Mx, delta_x, a);
end
DRP = U;


%% DRP-M格式
a = [-0.020843142770, 0.166705904415, -0.770882380518, 0, ...
     0.770882380518, -0.166705904415, 0.020843142770]; % j=-3到3
c_damp = [-0.014281184692, 0.086150669577, -0.235718815308, 0.327698660846, ...
          -0.235718815308, 0.086150669577, -0.014281184692]; % j=-3到3
nu_a = 0.001*delta_x^2;

% 使用RK4
for t = 1:3
    U1 = U-1/4*c*delta_t*(-1*DRPMcompute_rhs(U, Mx, delta_x, a, c, nu_a, c_damp));
    U2 = U-1/3*c*delta_t*(-1*DRPMcompute_rhs(U1, Mx, delta_x, a, c, nu_a, c_damp));
    U3 = U-1/2*c*delta_t*(-1*DRPMcompute_rhs(U2, Mx, delta_x, a, c, nu_a, c_damp));
    U = U-c*delta_t*(-1*DRPMcompute_rhs(U3, Mx, delta_x, a, c, nu_a, c_damp));
end
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
XXX = 0:delta_x:1;
plot(XXX,DRPM,'LineWidth',1.0);


%% 综合作图
exact = 1;
subplot(1,2,2);
XX = 0.8:delta_x:1;
X = 0.8:0.001:1;
for k = 1:64
    psik = r(k);
    Ek = (k/24)^4*exp(-2*(k/24)^2);
    exact = exact + 0.1*sqrt(Ek)*sin(2*pi*k*(X+psik));
end
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
exact = 1;
XX = 0:delta_x:1;
X = 0:0.001:1;
for k = 1:64
    psik = r(k);
    Ek = (k/24)^4*exp(-2*(k/24)^2);
    exact = exact + 0.1*sqrt(Ek)*sin(2*pi*k*(X+psik));
end
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

%% 误差计算
%精确解
X = 0:delta_x:1;
Exact = zeros(1,Mx+1);
for j = 1:Mx+1
    exact = 1;
    for k = 1:64
        psik = r(k);
        Ek = (k/24)^4*exp(-2*(k/24)^2);
        exact = exact + 0.1*sqrt(Ek)*sin(2*pi*k*(X(j)+psik));
    end
    Exact(j) = exact;
end

W1 = 0;
W2 = 0;
W3 = 0;
W4 = 0;
for j = 1:Mx
    %W1 = W1 + abs(DRP(j)-Exact(j))/Mx;
    W2 = W2 + abs(DRPM(j)-Exact(j))/Mx;
    %W3 = W3 + abs(MDCD(j)-Exact(j))/Mx;
    %W4 = W4 + abs(SADRP(j)-Exact(j))/Mx;
end

%% 误差作图
m = [0.1433 0.0647 0.0894 0.0895; 0.1019 0.0672 0.0761 0.0678; 0.0717 0.0380 0.0379 0.0179; 0.0113 0.0191 0.0070 6.65*10^(-4); 8.49*10^(-4) 0.0096 5.57*10^(-4) 1.64*10^(-5)];

y = [64 128 256 512 1024];

semilogy(y,m(:,1),'LineWidth',1.0);
hold on;
semilogy(y,m(:,2),'LineWidth',1.0);
semilogy(y,m(:,3),'LineWidth',1.0);
semilogy(y,m(:,4),'LineWidth',1.0);

xlabel('N');
ylabel('对数范数误差log(E)');
title('各格式的范数误差随空间网格数变化');
legend('DRP','DRP-M','MDCD','SA-DRP');


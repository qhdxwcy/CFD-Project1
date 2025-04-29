%% 4.3 初始条件
%clear;
Mx = 256;%空间网格数
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

%% 一阶迎风格式
n = 1;
% 使用RK4
for t = 1:Nt
    U1 = U-1/4*c*delta_t*UWcompute_dudx(U, Mx, delta_x, c, n);
    U2 = U-1/3*c*delta_t*UWcompute_dudx(U1, Mx, delta_x, c, n);
    U3 = U-1/2*c*delta_t*UWcompute_dudx(U2, Mx, delta_x, c, n);
    U = U-c*delta_t*UWcompute_dudx(U3, Mx, delta_x, c, n);
end
UW1 = U;

%% 二阶迎风格式
% 使用RK4
n = 2;
for t = 1:Nt
    U1 = U-1/4*c*delta_t*UWcompute_dudx(U, Mx, delta_x, c, n);
    U2 = U-1/3*c*delta_t*UWcompute_dudx(U1, Mx, delta_x, c, n);
    U3 = U-1/2*c*delta_t*UWcompute_dudx(U2, Mx, delta_x, c, n);
    U = U-c*delta_t*UWcompute_dudx(U3, Mx, delta_x, c, n);
end
UW2 = U;

%% 三阶迎风格式
% 使用RK4
n = 3;
for t = 1:Nt
    U1 = U-1/4*c*delta_t*UWcompute_dudx(U, Mx, delta_x, c, n);
    U2 = U-1/3*c*delta_t*UWcompute_dudx(U1, Mx, delta_x, c, n);
    U3 = U-1/2*c*delta_t*UWcompute_dudx(U2, Mx, delta_x, c, n);
    U = U-c*delta_t*UWcompute_dudx(U3, Mx, delta_x, c, n);
end
UW3 = U;

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
% exact = 0;
% for i = 1:20
%     exact = exact + sin(2*pi*i*(X-10));
% end
% exact = exact/20;
hold on;
plot(X,exact,'k');
scatter(XX,UW1(206:257),30,'bo');
scatter(XX,UW2(206:257),30,'rs');
scatter(XX,UW3(206:257),30,'gv');

plot(XX,UW1(206:257),'b','LineWidth',1.0);
plot(XX,UW2(206:257),'r','LineWidth',1.0);
plot(XX,UW3(206:257),'g','LineWidth',1.0);

xlabel('x','Fontsize',14);
ylabel('u','Fontsize',14);
legend('exact','一阶迎风','二阶迎风','三阶迎风','Location','best','Fontsize',14);
title('各阶迎风格式数值解（部分）','Fontsize',14);
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
% exact = 0;
% for i = 1:20
%     exact = exact + sin(2*pi*i*(X-10));
% end
% exact = exact/20;
hold on;
plot(X,exact,'k');
scatter(XX,UW1,10,'bo');
scatter(XX,UW2,10,'rs');
scatter(XX,UW3,10,'gv');

plot(XX,UW1,'b','LineWidth',1.0);
plot(XX,UW2,'r','LineWidth',1.0);
plot(XX,UW3,'g','LineWidth',1.0);

xlabel('x','Fontsize',14);
ylabel('u','Fontsize',14);
legend('exact','一阶迎风','二阶迎风','三阶迎风','Location','best','Fontsize',14);
title('各阶迎风格式数值解','Fontsize',14);
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
for j = 1:Mx
    W1 = W1 + abs(UW1(j)-Exact(j))/Mx;
    W2 = W2 + abs(UW2(j)-Exact(j))/Mx;
    %W3 = W3 + abs(UW3(j)-Exact(j))/Mx;
end
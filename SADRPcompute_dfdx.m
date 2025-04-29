%% SA-DRP格式空间离散计算函数
function dfdx = SADRPcompute_dfdx(f, N, dx, epsilon)
    dfdx = zeros(1,N+1);
    f_hat = zeros(1, N+1);
    for i = 1:N
        % 传感器
        f6 = zeros(1,6);
        for j = -2:3
            idx = mod(i +j - 1, N) + 1; % 周期边界处理
            f6(j+3) = f(idx);
        end
        S1 = f6(4)-2*f6(3)+f6(2);
        S2 = 1/4*(f6(5)-2*f6(3)+f6(1));
        S3 = f6(5)-2*f6(4)+f6(3);
        S4 = 1/4*(f6(6)-2*f6(4)+f6(2));
        C1 = f6(4)-f6(3);
        C2 = 1/3*(f6(5)-f6(2));
   
        k_ESW = acos(2*min((abs(abs(S1+S2)-abs(S1-S2))+abs(abs(S3+S4)-abs(S3-S4))+...
        abs(abs(C1+C2)-abs(C1-C2)/2)+2*epsilon)/(abs(S1+S2)+abs(S1-S2)+abs(S3+S4)+abs(S3-S4)+abs(C1+C2)+abs(C1-C2)+epsilon), 1)-1);

        if k_ESW >= 0 && k_ESW < 0.01
            a = 1/30;
        elseif k_ESW >= 0.01 && k_ESW < 2.5
            a = (k_ESW+1/6*sin(2*k_ESW)-4/3*sin(k_ESW))/(sin(3*k_ESW)-4*sin(2*k_ESW)+5*sin(k_ESW));
        else
            a = 0.1985842;
        end

        if k_ESW >= 0 && k_ESW <= 1
            b = 0.001;
        else
            b = 1*min(0.001+0.011*sqrt((k_ESW-1)/(pi-1)),0.012);
        end

        f_hat(i) = ((1/2*a+1/2*b)*f6(1)+(-3/2*a-5/2*b-1/12)*f6(2)+(a+5*b+7/12)*f6(3)+(a-5*b+7/12)*f6(4)+(-3/2*a+5/2*b-1/12)*f6(5)+(1/2*a-1/2*b)*f6(6));
        
    end
    f_hat(N+1) = f_hat(1);

    for i = 2:N
        dfdx(i) = (f_hat(i)-f_hat(i-1))/dx;
    end
    dfdx(1) = (f_hat(1)-f_hat(N))/dx;
    dfdx(N+1) = dfdx(1);

end

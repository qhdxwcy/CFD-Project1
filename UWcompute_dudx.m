%% 迎风格式空间离散计算函数
function dudx = UWcompute_dudx(U, N, dx, c, n)
    dudx = zeros(1, N+1);
    for i = 1:N
        u4 = zeros(1,4);
        for j = -2:1
            idx = mod(i +j - 1, N) + 1; % 周期边界处理
            u4(j+3) = U(idx);
        end
        if n == 1
            dudx(i) = c/dx*(u4(3)-u4(2));
        elseif n == 2
            dudx(i) = 1/2*c/dx*(3*u4(3)-4*u4(2)+u4(1));
        elseif n == 3
            dudx(i) = 1/6*c/dx*(2*u4(4)+3*u4(3)-6*u4(2)+u4(1));
        else
            disp('error');
        end
    end
    dudx(N+1) = dudx(1);
end
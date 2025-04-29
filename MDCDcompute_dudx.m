%% MDCD格式空间离散计算函数
function du_dx = MDCDcompute_dudx(U, N, dx, a, b)
    du_dx = zeros(1, N+1);
    s = zeros(1,7);
    s(1) = -1/2*a-1/2*b;
    s(2) = 2*a+3*b+1/12;
    s(3) = -5/2*a-15/2*b-2/3;
    s(4) = 10*b;
    s(5) = 5/2*a-15/2*b+2/3;
    s(6) = -2*a+3*b-1/12;
    s(7) = 1/2*a-1/2*b;
    for i = 1:N
        % 计算空间导数项
        du_dx(i) = 0;
        for j = -3:3
            idx = mod(i + j - 1, N) + 1; % 周期边界处理
            du_dx(i) = du_dx(i) + s(j+4) * U(idx);
        end
        du_dx(i) = du_dx(i) / dx;
    end
    du_dx(N+1) = du_dx(1);
end
%% DRP格式空间离散计算函数
function rhs = DRPcompute_rhs(U, N, dx, a)
    rhs = zeros(1, N+1);
    for i = 1:N
        % 计算空间导数项（DRP格式）
        du_dx = 0;
        for j = -3:3
            idx = mod(i + j - 1, N) + 1; % 周期边界处理
            du_dx = du_dx + a(j+4) * U(idx);
        end
        du_dx = du_dx / dx;
        
        rhs(i) = du_dx;
    end
    rhs(N+1) = rhs(1);
end
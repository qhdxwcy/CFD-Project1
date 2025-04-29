%% DRP-M格式右端项计算函数（含DRP格式和人工阻尼）
function rhs = DRPMcompute_rhs(U, N, dx, a, c, nu_a, c_damp)
    rhs = zeros(1, N+1);
    for i = 1:N
        % 计算空间导数项（DRP格式）
        du_dx = 0;
        for j = -3:3
            idx = mod(i + j - 1, N) + 1; % 周期边界处理
            du_dx = du_dx + a(j+4) * U(idx);
        end
        du_dx = du_dx / dx;
        
        % 计算人工阻尼项
        damp = 0;
        for j = -3:3
            idx = mod(i + j - 1, N) + 1;
            damp = damp + c_damp(j+4) * U(idx);
        end
        damp = -nu_a * damp / dx^2;
        
        % 合并右端项
        rhs(i) = -c * du_dx + damp;
    end
    rhs(N+1) = rhs(1);
end
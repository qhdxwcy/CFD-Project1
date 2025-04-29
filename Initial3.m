function u = Initial3(m,x)
    uu=0;
    for j = 1:m
        uu = uu+sin(2*pi*j*x);
    end
    u = uu/m;

end
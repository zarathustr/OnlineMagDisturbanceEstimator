function [c, ceq] = quat_con(x)
global mx my mz

mN = mN_f(x);
mD = mD_f(x);
k = x(4);
Mx = mx + x(1);
My = my + x(2);
Mz = mz + x(3);

c = [
    mN - 0.9;
    - mD - 0.8;
    1e-8 - abs(k);
    ];

ceq = Mx * Mx + My * My + Mz * Mz - 1;
end
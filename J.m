function f = J(x)
global ax ay az mx my mz q0 q1 q2 q3

mN = mN_f(x);
mD = mD_f(x);
k = x(4);
Mx = mx + x(1);
My = my + x(2);
Mz = mz + x(3);

f = (ax*My - ay*Mx - ay*mN - k*q0).^2 + ...
         (-Mx + az*Mx - ax*Mz + (az - 1)*mN + ax*mD - k*q1).^2 + ...
         (-My + az*My - ay*Mz + ay*mD - k*q2).^2 + ...
         (-Mz + az*mD - ax*mN - k*q3).^2;
end
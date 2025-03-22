function [A, B] = matrices(x_ss, u_ss, p)

syms x [p.n_x 1]
syms u [p.n_u 1]
dxdt = dynamics(x, u, p);
dfdx = matlabFunction(jacobian(dxdt, x), 'Vars', {x, u});
dfdu = matlabFunction(jacobian(dxdt, u), 'Vars', {x, u});

A = cell(p.n_equi, 1); 
B = cell(p.n_equi, 1); 

for i = 1:p.n_equi
    A{i} = dfdx(x_ss{i}, u_ss{i});
    B{i} = dfdu(x_ss{i}, u_ss{i});
    [A{i}, B{i}] = c2d(A{i}, B{i}, p.ts);
end

end
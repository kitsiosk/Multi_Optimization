syms f(x, y) g(x, y);
f(x, y) = (x^3) * exp(-x^2 - y^4);
g(x, y) = x^4 + y^2 - 0.2*sin(2*pi*x) - 0.3*cos(2*pi*y);

% Current function. Switch between f and g
syms t(x, y);
t(x, y) = f(x, y);

% 3) Newton method
epsilon = 0.1;     % Termination constant
k = 0;               % Number of repetitions
gamma0 = 0.1;       % constant | minimizing f(x + g*d) | Armijo rule
% Booleans to define the choice of gamma
isConstant = true;
isMinimizing = false;
isArmijo = false;

% Initialize x0, y0
xk = -1; % 0 | -1 | 1
yk = -1; % 0 | -1 | 1

grad = gradient(t);
hess = hessian(t);

while norm(grad(xk, yk)) > epsilon
    % Use vpa() to convert sin(), cos() and  pi tofloating point approximations
    % with 10 decimals in order to spped up calculations
    hessk = hess(xk, yk);
    gradk = grad(xk, yk);
    eigenvalues = eig(hessk);
    if ~ (all(eigenvalues > 0) || all(eigenvalues < 0))
        fprintf('Hessian not positive(negative) definite. Newton method cannot apply\n');
        return;
    end
    if det(hessk) == 0
        fprintf('Hessian not invertible, aborting\n');
        return;
    end
    
    if isConstant
        gammak = gamma0;
    elseif isMinimizing
        syms p(u)
        p(u) = f(xk + u*dk(1), yk + u*dk(2));
        grad_p = vpa(diff(p, u));
        a = 0;
        b = 0.1;
        l = 0.001;
        n = ceil(-log2( l/(b-a)));

        for i = 1:n
            x0 = (a + b)/2;
            if grad_p(x0) == 0
                a=x0;
                b=x0;
                break;
            elseif grad_p(x0) > 0
                a = a;
                b = x0;
            else
                a = x0;
                b = b;
            end
        end
        gammak = x0;
    elseif isArmijo
        m=0;
        gammak = gamma0;
        xk_next = xk + gamma0 * dk(1);
        yk_next = yk + gamma0 * dk(2);
        while f(xk, yk) - f(xk_next, yk_next) < a*exp(-m)*gamma0*(norm(dk)^2)
            gammak = gamma0*exp(-m);
            xk_next = xk + gammak * dk(1);
            yk_next = yk + gammak * dk(2);
            m = m+1;
        end
    end
    
    dk = hessk\gradk;
    dk = -gammak*dk; % direction vector
    xk = xk + dk(1);
    yk = yk + dk(2);
    k = k + 1;
end

fprintf("Potential local optimum found at (%f, %f)\n", xk, yk);
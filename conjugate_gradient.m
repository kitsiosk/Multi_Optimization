syms f(x, y) g(x, y);
f(x, y) = (x^3) * exp(-x^2 - y^4);
g(x, y) = x^4 + y^2 - 0.2*sin(2*pi*x) - 0.3*cos(2*pi*y);

% Current function. Switch between f and g
syms t(x, y);
t(x, y) = g(x, y);

% 5) Conjugate gradient method
epsilon = 0.01;     % Termination constant
k = 0;             % Number of repetitions
gamma0 = 0.01;

% Booleans to define the choice of gamma
isConstant = false;
isMinimizing = true;

% Initialize x0, y0
xk = -0.6; % 0 | -0.6 | 1
yk = -0.6; % 0 | -0.6 | 1
% Initialize X(k-1) that holds the previous value of xk
xk_prev = xk;
yk_prev = yk;

values = [g(xk, yk)];

grad = gradient(t);

while norm(grad(xk, yk)) > epsilon
    % Use vpa() to convert sin(), cos() and  pi tofloating point approximations
    % with 10 decimals in order to speed up calculations
    if k==0
        dk = -vpa(grad(xk, yk), 10);
        dk_prev = dk;
    else
        bk = ((grad(xk, yk)')*(grad(xk, yk) - grad(xk_prev, yk_prev)))/((grad(xk_prev, yk_prev)')*grad(xk_prev, yk_prev));
        dk = vpa(-grad(xk, yk) + bk*dk_prev);
    end
    
    if isConstant
        gammak = gamma0;
    elseif isMinimizing
        syms p(u)
        p(u) = g(xk + u*dk(1), yk + u*dk(2));
        grad_p = vpa(diff(p, u));
        a = 0;
        b = 0.5;
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
    end
    
    xk_prev = xk;
    yk_prev = yk;
    xk = xk + gammak * dk(1);
    yk = yk + gammak * dk(2);
    k = k + 1;
    values = [values; g(xk, yk)];
end

fprintf("Potential local optimum found at (%f, %f)\n", xk, yk);
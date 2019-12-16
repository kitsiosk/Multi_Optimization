syms f(x, y) g(x, y);
f(x, y) = (x^3) * exp(-x^2 - y^4);
g(x, y) = x^4 + y^2 - 0.2*sin(2*pi*x) - 0.3*cos(2*pi*y);

% Current function. Switch between f and g
syms t(x, y);
t(x, y) = f(x, y);

% 1) Plot the functions in 3D space
plot = false;
if plot
    figure(1);
    fsurf(f);
    title('f(x, y)');
    xlabel('x');
    ylabel('y');
    zlabel('z');

    figure(2);
    fsurf(g);
    title('g(x, y), xe[-5, 5], ye[-5, 5]');
    xlabel('x');
    ylabel('y');
    zlabel('z');

    figure(3);
    fsurf(g, [-1 1 -1 1]);
    title('g(x, y), xe[-1, 1], ye[-1, 1]');
    xlabel('x');
    ylabel('y');
    zlabel('z');
end

% 2) Steepest descent algorithm
epsilon = 0.01;     % Termination constant
k = 0;               % Number of repetitions
gamma0 = 0.01;       % constant | minimizing f(x + g*d) | Armijo rule
% Booleans to define the choice of gamma
isConstant = false;
isMinimizing = false;
isArmijo = true;

% Initialize x0, y0
xk = 1; % 0 | -1 | 1
yk = 1; % 0 | -1 | 1

values = [f(xk, yk)];

grad = gradient(t);

while norm(grad(xk, yk)) > epsilon
    % Use vpa() to convert sin(), cos() and  pi tofloating point approximations
    % with 10 decimals in order to spped up calculations
    dk = -vpa(grad(xk, yk), 10);
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
    
    xk = xk + gammak * dk(1);
    yk = yk + gammak * dk(2);
    k = k + 1;
    values = [values; f(xk, yk)];
end

fprintf("Potential local optimum found at (%f, %f)\n", xk, yk);
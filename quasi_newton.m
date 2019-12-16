syms f(x, y) g(x, y);
f(x, y) = (x^3) * exp(-x^2 - y^4);
g(x, y) = x^4 + y^2 - 0.2*sin(2*pi*x) - 0.3*cos(2*pi*y);

% Current function. Switch between f and g
syms t(x, y);
t(x, y) = g(x, y);

% 3) Newton method
epsilon = 0.01;     % Termination constant
k = 0;               % Number of repetitions
gamma0 = 0.1;
ksi_k = 1;          % Value of î
ee = 0;      % an epsilon value to avoid division by 0

% Booleans to define the choice of gamma
isConstant = true;
isMinimizing = true;

% Initialize x0, y0 and D0
xk = 0; % 0 | -0.6 | 1
yk = 0; % 0 | -0.6 | 1
xk_prev = xk;
yk_prev = yk;
D0 = eye(2);

values = [];

grad = gradient(t);

while norm(grad(xk, yk)) > epsilon
    % Use vpa() to convert sin(), cos() and  pi to floating point approximations
    % with 10 decimals in order to spped up calculations
    gradk = vpa(grad(xk, yk));
    
    if k==0
        Dk = D0;
    else
        pk = xk_prev - xk;
        qk = vpa(grad(xk, yk) - grad(xk_prev, yk_prev));
        tk = (qk')*Dk*qk;
        uk = pk/((pk') * qk + ee) - (Dk*qk)/(tk + ee);
        Dk = Dk + ((pk')*pk)/((pk')*qk + ee) - (Dk*qk*(qk')*Dk)/((qk')*Dk*qk + ee) + ...
            + ksi_k*tk*uk*(uk');
    end

    if isConstant
        gammak = gamma0;
    elseif isMinimizing
        syms p(u)
        p(u) = t(xk + u*dk(1), yk + u*dk(2));
        grad_p = vpa(diff(p, u));
        a = 0;
        b = 0.5;
        l = 0.001;
        n = ceil(-log2( l/(b-a)));

        for k = 1:n
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
    
    values = [values; g(xk, yk)];

    dk = -Dk*gradk;
    xk = xk + gammak*dk(1);
    yk = yk + gammak*dk(2);
    k = k + 1;
    
    % exit if k>n
    if k>2
        break;
    end
end

fprintf("Potential local optimum found at (%f, %f)\n", xk, yk);
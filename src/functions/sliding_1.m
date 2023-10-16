function [] = sliding_1(x0)

    clc

    k = 5;
    m = 1;
    u = 0.75;
    g = 9.81;
    w = sqrt(k/m);

    tspan = [0 15];
    opts = odeset('Refine', 10);
    %opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
    sol = ode15s(@odefun, tspan, x0, opts);
    % plot(sol.x, sol.y(2,:))

    t = linspace(0, 15 , 1000);
    x = deval(sol, t, 1);
    y = deval(sol, t, 2);

    % Position and Velocity in time
    figure()
    plot(t, x);
    hold on
    xd =  @(t) 2+sin(2.5*t)+2*cos(1.25*t);
    plot(t, xd(t));
    title('Sliding1: Trajectory Tracking')
    xlabel('{t} (sec)')
    ylabel('{x_1} (m)')
    legend('System Trajectory ', 'Desired Trajectory')
    hold off

    function dxdt = odefun(t,x)

        k = 5;
        m = 1;
        u = 0.75;
        g = 9.81;
        w = sqrt(k/m);

        m_approx = 1.25;
        u_approx = 0.625;
        k_approx = 7;
        Dk_approx = 0.625;
        a_approx = 0.75;
        b_approx = 0.03;

        lamda = 15;
        c = 100;

        f = @(x, t) (5+0.6*exp(-0.045*t))*x*(1+0.9*x^2);

        sat = @(x, delta) min(max(x/delta, -1), 1);

        xd = @(t) 2+sin(2.5*t)+2*cos(1.25*t);
        xdd = @(t) 2.5*cos(2.5*t)-2.5*sin(1.25*t);
        xddd = @(t) -6.25*sin(2.5*t)-3.125*cos(1.25*t);

        err = x(1) - xd(t);
        errd = w*x(2) - xdd(t);

        s = errd + lamda*err;

        dxdt = zeros(2,1);
        dxdt(1) = w*x(2);

        r = 0.9375*g*abs(sat(w*x(2), 1)) + 0.75*abs(xddd(t)-lamda*errd) + 3.375*abs(x(1)) + 4.4375*(abs(x(1)))^3 + c;
        control = u_approx*m_approx*g*sat(w*x(2), 1) + m_approx*(xddd(t)-lamda*errd) - r*sat(s, 1) ...
            + (k_approx + Dk_approx*(1-exp(-b_approx*t)))*x(1)*(1+a_approx*(x(1))^2);

        if x(2)>1e-8
            dxdt(2) = -u*g/w-1/(w*m)*f(x(1),t)+1/(w*m)*control;
        elseif x(2)<-1e-8
            dxdt(2) = +u*g/w-1/(w*m)*f(x(1),t)+1/(w*m)*control;
        else
            if f(x(1),t)<-u*m*g + control
                dxdt(2) = -u*g/w-1/(w*m)*f(x(1),t)+1/(w*m)*control;
            elseif f(x(1),t)>u*m*g + control
                dxdt(2) = +u*g/w-1/(w*m)*f(x(1),t)+1/(w*m)*control;
            else
                dxdt(2) = 0;
            end
        end

    end

end


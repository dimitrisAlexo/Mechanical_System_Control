function [] = redesign_2(x0)

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
    title('Redesign2: Trajectory Tracking')
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

        a1 = -50;  % eigenvalues
        a2 = -50;

        l1 = m*a1*a2;
        l2 = -m*(a1+a2);
        F = [l1; l2];
        b = [0; 1];

        p11 = l2/(2*l1) + (l1*m+l1^2)/(2*l1*l2);
        p12 = m/(2*l1);
        p22 = (m^2+m*l1)/(2*l1*l2);
        P = [p11 p12; p12 p22];

        f = @(x, t) (5+0.6*exp(-0.045*t))*x*(1+0.9*x^2);

        sat = @(x, delta) min(max(x/delta, -1), 1);

        xd = @(t) 2+sin(2.5*t)+2*cos(1.25*t);
        xdd = @(t) 2.5*cos(2.5*t)-2.5*sin(1.25*t);
        xddd = @(t) -6.25*sin(2.5*t)-3.125*cos(1.25*t);

        dxdt = zeros(2,1);
        dxdt(1) = w*x(2);
        
        e = [x(1) - xd(t); dxdt(1) - xdd(t)];

        control = u_approx*m_approx*g*sat(w*x(2), 1) + m*xddd(t) ...
            - transpose(F)*e - (0.9375*g*abs(sat(w*x(2), 1)) + 0.75*abs(xddd(t)) ...
            + 22*(abs(x(1)))^3)*sat(transpose(b)*P*e, 1);

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


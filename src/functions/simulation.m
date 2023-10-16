function [] = simulation(x0)

    clc

    k = 5;
    m = 1;
    u = 0.75;
    g = 9.81;
    w = sqrt(k/m);

    tspan = [0 5];
    opts = odeset('Refine', 10);
    sol = ode15s(@odefun, tspan, x0, opts);

    t = linspace(0, 5 , 2000);
    x = deval(sol, t, 1);
    y = deval(sol, t, 2);

    % Position and Velocity in time
    figure()
    plot(t, x);
    hold on
    plot(t, y);
    title('State Variables')
    xlabel('{t} (sec)')
    ylabel('{x_1, x_2}')
    legend('Position (m)', 'Velocity (m/s)')
    hold off

    % Phase Depiction
    figure()
    hold on
    axis square
    title('Phase Depiction')
    xlabel('{x_1} (m)')
    ylabel('{x_2} (m/s)')
    xlim([-5 5])
    ylim([-5 5])
    plot(x, y);
    plot([-5 5], [0 0], 'k-', [0 0], [-5 5], 'k', 'linewidth', 1)
    hold off

    for t = 0:0.01:5
        x = deval(sol, t, 1);
        y = deval(sol, t, 2);
        if (x < u*m*g/k && x > -u*m*g/k && y < 1e-8 && y > -1e-8)
            disp("Resting time = " + t + "sec");
            break;
        end
    end

    function dxdt = odefun(~,x)

        k = 5;
        m = 1;
        u = 0.75;
        g = 9.81;
        w = sqrt(k/m);

        dxdt = zeros(2,1);
        dxdt(1) = w*x(2);

        if x(2)>1e-8
            dxdt(2) = -u*g/w-k/(w*m)*x(1);
        elseif x(2)<-1e-8
            dxdt(2) = +u*g/w-k/(w*m)*x(1);
        else
            if x(1)<-u*m*g/k
                dxdt(2) = -u*g/w-k/(w*m)*x(1);
            elseif x(1)>u*m*g/k
                dxdt(2) = +u*g/w-k/(w*m)*x(1);
            else
                dxdt(2) = 0;
            end
        end

    end

end

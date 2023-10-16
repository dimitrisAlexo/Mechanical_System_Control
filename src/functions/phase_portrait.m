function [] = phase_portrait()

    clc

    k = 5;
    m = 1;
    u = 0.75;
    g = 9.81;
    w = sqrt(k/m);

    figure();
    axis square
    hold on
    title('Phase Portrait')
    xlabel('{x_1} (m)')
    ylabel('{x_2} (m/s)')

    for c = [-6.7 0 10 23]
        fimplicit(@(x_1, x_2) k/m*x_1.^2 + w^2*x_2.^2+2*u*g*x_1-c, [-5 5 0 5])
        fimplicit(@(x_1, x_2) k/m*x_1.^2 + w^2*x_2.^2-2*u*g*x_1-c, [-5 5 -5 0])
    end

    plot([-5 5], [0 0], 'k-', [0 0], [-5 5], 'k', 'linewidth', 1)
    hold off

end


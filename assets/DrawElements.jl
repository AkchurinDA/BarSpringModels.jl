using Luxor
using MathTeXEngine

# Connection element:
@drawsvg begin
    # Background:
    @layer begin
        setcolor("white")
        box(O, 600, 600, 50, :fill)
        setcolor("black")
        box(O, 599, 599, 47.5, :stroke)
    end

    # Connection element:
    @layer begin
        # Verticals:
        line(Point(-100, +200), Point(-100, -200), action=:stroke)
        line(Point(+100, +200), Point(+100, -200), action=:stroke)

        # Tranlational springs:
        line(Point(-100, -150), Point(-50, -150), action=:stroke)
        line(Point(+50, -150), Point(+100, -150), action=:stroke)
        Points = [Point(x, -150 + 10 * sin(5 * x * pi / 50)) for x in -50:0.01:+50]
        poly(Points, action=:stroke)

        line(Point(-100, -100), Point(-50, -100), action=:stroke)
        line(Point(+50, -100), Point(+100, -100), action=:stroke)
        Points = [Point(x, -100 + 10 * sin(5 * x * pi / 50)) for x in -50:0.01:+50]
        poly(Points, action=:stroke)

        line(Point(-100, -50), Point(-50, -50), action=:stroke)
        line(Point(+50, -50), Point(+100, -50), action=:stroke)
        Points = [Point(x, -50 + 10 * sin(5 * x * pi / 50)) for x in -50:0.01:+50]
        poly(Points, action=:stroke)

        # Rotational springs:
        line(Point(-100, +50), Point(0, +50), action=:stroke)
        line(Point(+18, +50), Point(+100, +50), action=:stroke)
        gsave()
        translate(0, +50)
        spiral(1, 1, action=:stroke, period=6 * pi)
        grestore()

        line(Point(-100, +100), Point(0, +100), action=:stroke)
        line(Point(+18, +100), Point(+100, +100), action=:stroke)
        gsave()
        translate(0, +100)
        spiral(1, 1, action=:stroke, period=6 * pi)
        grestore()

        line(Point(-100, +150), Point(0, +150), action=:stroke)
        line(Point(+18, +150), Point(+100, +150), action=:stroke)
        gsave()
        translate(0, +150)
        spiral(1, 1, action=:stroke, period=6 * pi)
        grestore()
    end

    # DOFs:
    @layer begin
        circle(Point(-100, 0), 5, action=:fill)
        circle(Point(+100, 0), 5, action=:fill)
    end

    # Labels:
    @layer begin
        fontsize(15)
        text(L"$i$", Point(-100 - 10, 0), halign=:right, valign=:middle)
        text(L"$j$", Point(+100 + 10, 0), halign=:left, valign=:middle)

        text(L"$u_{x_i}$", Point(-100 + 10, -150 - 10), halign=:left, valign=:bottom)
        text(L"$u_{x_j}$", Point(+100 - 10, -150 - 10), halign=:right, valign=:bottom)

        text(L"$u_{y_i}$", Point(-100 + 10, -100 - 10), halign=:left, valign=:bottom)
        text(L"$u_{y_j}$", Point(+100 - 10, -100 - 10), halign=:right, valign=:bottom)

        text(L"$u_{z_i}$", Point(-100 + 10, -50 - 10), halign=:left, valign=:bottom)
        text(L"$u_{z_j}$", Point(+100 - 10, -50 - 10), halign=:right, valign=:bottom)

        text(L"$\theta_{x_i}$", Point(-100 + 10, +50 - 10), halign=:left, valign=:bottom)
        text(L"$\theta_{x_j}$", Point(+100 - 10, +50 - 10), halign=:right, valign=:bottom)

        text(L"$\theta_{y_i}$", Point(-100 + 10, +100 - 10), halign=:left, valign=:bottom)
        text(L"$\theta_{y_j}$", Point(+100 - 10, +100 - 10), halign=:right, valign=:bottom)

        text(L"$\theta_{z_i}$", Point(-100 + 10, +150 - 10), halign=:left, valign=:bottom)
        text(L"$\theta_{z_j}$", Point(+100 - 10, +150 - 10), halign=:right, valign=:bottom)
    end
end

# Rigid link:
@drawsvg begin
    # Background:
    @layer begin
        setcolor("white")
        box(O, 600, 600, 50, :fill)
        setcolor("black")
        box(O, 599, 599, 47.5, :stroke)
    end

    # Axes:
    @layer begin
        arrow(Point(-150, +150), Point(+200, +150))
        arrow(Point(-150, +150), Point(-150, -200))
        arrow(Point(-150, +150), Point(-250, +250))
    end

    # Rigid link:
    @layer begin
        line(Point(-100, +100), Point(+150, 0), action=:stroke)
    end

    # DOFs:
    @layer begin
        # Master node:
        arrow(Point(-100, +100), Point(-100 + 50, +100))
        arrow(Point(-100, +100), Point(-100, +100 - 50))
        arrow(Point(-100, +100), Point(-100 - 100 * sqrt(2) / 7, +100 + 100 * sqrt(2) / 7))

        # Slave node:
        arrow(Point(+150, 0), Point(+150 + 50, 0))
        arrow(Point(+150, 0), Point(+150, 0 - 50))
        arrow(Point(+150, 0), Point(+150 - 100 * sqrt(2) / 7, 0 + 100 * sqrt(2) / 7))
    end

    # Labels:
    @layer begin
        fontsize(15)
        text(L"$x$", Point(+200 + 5, +150), halign=:left, valign=:middle)
        text(L"$y$", Point(-150, -200 - 5), halign=:center, valign=:bottom)
        text(L"$z$", Point(-250 - 5, +250 + 5), halign=:right, valign=:top)
        text("Master node", Point(-100, +100 + 10), halign=:left, valign=:top)
        text(L"$u_{x_M}, \theta_{x_M}$", Point(-100 + 50 + 5, +100), halign=:left, valign=:middle)
        text(L"$u_{y_M}, \theta_{y_M}$", Point(-100, +100 - 50 - 5), halign=:center, valign=:bottom)
        text(L"$u_{z_M}, \theta_{z_M}$", Point(-100 - 100 * sqrt(2) / 7, +100 + 100 * sqrt(2) / 7 + 5), halign=:center, valign=:top)
        text("Slave node", Point(+150, 0 + 10), halign=:left, valign=:top)
        text(L"$u_{x_S}, \theta_{x_S}$", Point(+150 + 50 + 5, 0), halign=:left, valign=:middle)
        text(L"$u_{y_S}, \theta_{y_S}$", Point(+150, 0 - 50 - 5), halign=:center, valign=:bottom)
        text(L"$u_{z_S}, \theta_{z_S}$", Point(+150 - 100 * sqrt(2) / 7, 0 + 100 * sqrt(2) / 7 + 5), halign=:center, valign=:top)
    end
end
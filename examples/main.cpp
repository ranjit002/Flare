#include <SFML/Graphics.hpp>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>

#include "Flare/boundary.h"
#include "Flare/fluid.h"
#include "Flare/solver.h"

// Convert velocity & density to color
sf::Color colorFromVelDensity(const Eigen::Vector3f& vel,
    float maxSpeed,
    float density,
    float maxDensity);

int main()
{
    // Fluid dimensions
    const int W = 64, H = 64, D = 64;
    const float dt = 0.1f;
    const float maxDensity = 100.f;
    const float maxSpeed = 30.f;

    // Target window size
    const int targetWindowSize = 1000;
    const int scaleX = targetWindowSize / W;
    const int scaleY = targetWindowSize / H;
    const int SCALE_FACTOR = std::max(1, std::min(scaleX, scaleY));

    fluid::Fluid fluid{W, H, D, 10};
    solver::BasicSolver solver{10.0f, 1.f, 10, 10};

    // Boundaries
    std::vector<std::unique_ptr<boundary::IBoundary>> bcs;
    bcs.push_back(std::make_unique<boundary::Inflow>(
        Eigen::Vector3f(maxSpeed, 0, 0), 100));
    bcs.push_back(std::make_unique<boundary::Circle>(
        Eigen::Vector3f(W / 2.f, H / 2.f, D / 2.f), 10));
    bcs.push_back(std::make_unique<boundary::Box>(W, H, D));

    for (auto& bc : bcs) solver.addBC(std::move(bc));

    // SFML window
    sf::RenderWindow window(sf::VideoMode({SCALE_FACTOR * W, SCALE_FACTOR * H}),
        "Fluid Simulation");
    window.setFramerateLimit(60);

    sf::RectangleShape pixel(sf::Vector2f(SCALE_FACTOR, SCALE_FACTOR));

    while (window.isOpen())
    {
        while (const std::optional<sf::Event> event = window.pollEvent())
            if (event->is<sf::Event::Closed>())
                window.close();

        solver.step(fluid, dt);
        window.clear();

        for (int y = 0; y < H; ++y)
        {
            for (int x = 0; x < W; ++x)
            {
                bool isBorder = (x == 0 || x == W - 1 || y == 0 || y == H - 1);
                bool isSolid = false;
                if (!isBorder)
                {
                    for (auto& bc : solver.boundaries())
                        if (bc->isSolid(x, y, D / 2))
                        {
                            isSolid = true;
                            break;
                        }
                }

                sf::Color fillColor =
                    isSolid
                        ? sf::Color::White
                        : colorFromVelDensity(fluid.getVelocity(x, y, D / 2),
                              maxSpeed,
                              fluid.getDensity(x, y, D / 2),
                              maxDensity);

                pixel.setFillColor(fillColor);
                pixel.setPosition(
                    {x * float(SCALE_FACTOR), y * float(SCALE_FACTOR)});
                window.draw(pixel);
            }
        }

        window.display();
    }
}

sf::Color colorFromVelDensity(const Eigen::Vector3f& vel,
    float maxSpeed,
    float density,
    float maxDensity)
{
    float speed = vel.norm();
    float t = std::clamp(speed / maxSpeed, 0.f, 1.f);
    float dens = std::clamp(density / maxDensity, 0.f, 1.f);

    float hue = (1.f - t) * 240.f;
    float c = 1.f;
    float h = hue / 60.f;
    float xC = c * (1 - fabs(fmod(h, 2.f) - 1.f));

    float rF, gF, bF;
    if (h < 1)
    {
        rF = c;
        gF = xC;
        bF = 0;
    }
    else if (h < 2)
    {
        rF = xC;
        gF = c;
        bF = 0;
    }
    else if (h < 3)
    {
        rF = 0;
        gF = c;
        bF = xC;
    }
    else if (h < 4)
    {
        rF = 0;
        gF = xC;
        bF = c;
    }
    else if (h < 5)
    {
        rF = xC;
        gF = 0;
        bF = c;
    }
    else
    {
        rF = c;
        gF = 0;
        bF = xC;
    }

    return sf::Color(static_cast<std::uint8_t>(rF * 255.f * dens),
        static_cast<std::uint8_t>(gF * 255.f * dens),
        static_cast<std::uint8_t>(bF * 255.f * dens));
}
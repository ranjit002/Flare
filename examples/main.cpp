#include <SFML/Graphics.hpp>

#include "Flare/fluid/Fluid.h"
#include "Flare/solver/BasicSolver.h"

int main()
{
  int FRAME_RATE = 1;
  const int W = 32, H = 32, D = 2;
  const float dt = 0.1f;

  fluid::Fluid fluid(W, H, D);
  fluid.setDensity(W / 2, H / 2, D / 2, 100.0f);
  solver::BasicSolver solver;

  sf::RenderWindow window(sf::VideoMode({W, H}), "Fluid Simulation");
  window.setFramerateLimit(FRAME_RATE);

  while (window.isOpen())
  {
    while (const std::optional<sf::Event> event = window.pollEvent())
      if (event->is<sf::Event::Closed>()) window.close();

    window.clear();
    solver.step(fluid, dt);

    float maxDensity = 100.f;
    sf::RectangleShape pixel(sf::Vector2f(1.f, 1.f));

    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x)
      {
        float val = fluid.getDensity(x, y, D / 2) / maxDensity;
        pixel.setFillColor(sf::Color(static_cast<int>(val * 255.f),
            static_cast<int>(val * 255.f),
            static_cast<int>(val * 255.f)));
        pixel.setPosition({x * pixel.getSize().x, y * pixel.getSize().y});
        window.draw(pixel);
      }

    window.display();
  }
}
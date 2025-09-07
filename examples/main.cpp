#include <SFML/Graphics.hpp>

#include "Flare/fluid/Fluid.h"
#include "Flare/solver/BasicSolver.h"

int main()
{
  const int SCALE_FACTOR = 10;
  const int FRAME_RATE = 30;
  const int W = 32, H = 32, D = 2;
  const float dt = 0.1f;

  fluid::Fluid fluid{W, H, D};
  fluid.setDensity(W / 2, H / 2, D / 2, 100.0f);
  solver::BasicSolver solver{};

  sf::RenderWindow window(
      sf::VideoMode({SCALE_FACTOR * W, SCALE_FACTOR * H}), "Fluid Simulation");
  window.setFramerateLimit(FRAME_RATE);

  while (window.isOpen())
  {
    while (const std::optional<sf::Event> event = window.pollEvent())
      if (event->is<sf::Event::Closed>()) window.close();

    window.clear();
    solver.step(fluid, dt);

    float maxDensity = 100.f;
    sf::RectangleShape pixel(sf::Vector2f(SCALE_FACTOR, SCALE_FACTOR));

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
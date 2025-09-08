#include <SFML/Graphics.hpp>
#include <memory>

#include "Flare/boundary/Box.h"
#include "Flare/boundary/Circle.h"
#include "Flare/boundary/IBoundary.h"
#include "Flare/boundary/Inflow.h"
#include "Flare/fluid/Fluid.h"
#include "Flare/solver/BasicSolver.h"

void addDensityInflow(fluid::Fluid& fluid, int H, int D)
{
  for (int y = 0; y < H; ++y)
    for (int z = 0; z < D; ++z) fluid.setDensity(0, y, z, 100.0f);
}

int main()
{
  const int SCALE_FACTOR = 10;
  const int FRAME_RATE = 5;
  const int W = 32, H = 32, D = 32;
  const float dt = 0.1f;

  fluid::Fluid fluid{W, H, D};
  solver::BasicSolver solver{0.0f, 10.f, 10, 10};

  std::vector<std::unique_ptr<boundary::IBoundary>> bcs;

  Eigen::Vector3f inflowVel(100.0f, 100.0f, 100.0f);
  bcs.push_back(std::make_unique<boundary::InflowBoundary>(inflowVel));

  bcs.push_back(std::make_unique<boundary::BoxBoundary>(W, H, D));

  Eigen::Vector3f centre = Eigen::Vector3f{static_cast<float>(W) / 2,
      static_cast<float>(H) / 2,
      static_cast<float>(D) / 2};

  bcs.push_back(std::make_unique<boundary::CircleBoundary>(centre, 5));

  for (auto& bc : bcs)
  {
    solver.addBC(std::move(bc));
  }

  sf::RenderWindow window(
      sf::VideoMode({SCALE_FACTOR * W, SCALE_FACTOR * H}), "Fluid Simulation");
  window.setFramerateLimit(FRAME_RATE);

  while (window.isOpen())
  {
    while (const std::optional<sf::Event> event = window.pollEvent())
      if (event->is<sf::Event::Closed>()) window.close();

    window.clear();
    addDensityInflow(fluid, H, D);
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
# Flare

<table>
<tr>
<td width="160" valign="middle">
  <img src="./assets/flare.png" alt="Flare Logo" width="150"/>
</td>
<td valign="middle">
  <p>
    <b>Flare</b> is a <b>high-performance 3D fluid simulation library</b> built in C++.  
    It focuses on <b>speed, efficiency, and practical simulation features</b> while remaining modular enough for experimentation and extension.
  </p>
</td>
</tr>
</table>

---

## Features

* **Structure of Arrays (SoA)**: Optimizes memory access patterns and improves SIMD-friendliness for velocity and density fields.
* **Multithreading**: Uses OpenMP to parallelize solver loops, taking full advantage of multicore CPUs.
* **Semi-Lagrangian Advection**: Stable method for transporting velocity and density fields even at larger time steps.
* **Customizable Boundaries**: Predefined walls, inflows, and obstacles designed for minimal overhead.

---

## Getting Started

### Prerequisites

* C++17 or newer
* [CMake 3.15+](https://cmake.org/)
* [Eigen](https://eigen.tuxfamily.org/)
* SFML
* OpenMP


## Quick Example

```cpp
#include <Flare/fluid.h>
#include <Flare/solver.h>
#include <Flare/boundary.h>

int main() {
    fluid::Fluid fluid(32, 32, 32);
    solver::BasicSolver solver;

    solver.addBC(boundary::Box(32, 32, 32));

    const float dt = 0.1f;
    for (int i = 0; i < 100; ++i) 
    {
        solver.step(fluid, dt);
    }
    return 0;
}
```

---

## Roadmap

* GPU acceleration (CUDA / OpenCL)
* More advanced solver (e.g. Finite element methods)
* More boundary types and inflow/outflow effects
* Built-in visualization tools for debugging and demoing fluid behavior

---
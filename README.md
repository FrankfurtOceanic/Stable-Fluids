# Stable-Fluids
 
This is a C++ implementation of Stam's Stable fluids. The paper's methods were extended to include a buoyancy step. 
* Functionality: a 2d Eulerian fluid based on Stable Fluids (2003), Stam et al.
* Programmed methods for adding density sources, diffusion, advection and projections to maintain the incompressibility constraint. Additionally, buoyancy we added to emulate the effect of warm air rising.
* The system works on a scalar level to move densities/temperatures within the system as well as modeling the effects on velocity.


Here's a video showing off the end results where the velocities of each grid cell is shown:
[![Fluids](https://img.youtube.com/vi/a83um4R1h2A/0.jpg)](https://youtu.be/a83um4R1h2A "Fluids")  

https://youtu.be/a83um4R1h2A



Please use CMake to build the project for maximum compatibility. Some third party libraries are required. They are eigen, freetype, glew, glw, glm, imgui, and rapidxml.


# TrivialConnections
Discrete Differential Geometry Exterior Calculus:  Keenan Crane's Trivial Connections Paper

### Some modifications(Make this repository buildable on Windows)

This fork has made some modifications based on the original repository and can now be successfully configured on Windows 11. 

Please update the `SuiteSparse` configuration path in CMakeLists.txt and build using CMake.

### Usage, e.g.

make <br/>
./connection ./meshes/torus.obj

### converting the code from mac to linux:

 suitesparse Changed from this: <br/>
   #include <SuiteSparseQR.hpp> <br/>
   #include <umfpack.h> <br/>
   #include <cholmod.h>

 to this: <br/>
   #include <suitesparse/SuiteSparseQR.hpp> <br/>
   #include <suitesparse/umfpack.h>
   
Change UF_long to SuiteSparse_long <br/>
   #include <suitesparse/cholmod.h>

lso, 
 I needed to add "-lumfpack -lamd"
 to the linker for umf libs.

 GL changed from this<br/>
   #include <GLUT/glut.h>

to this  <br/>
   #include <GL/glut.h>

And added glew include and "-lGLEW " link  <br/>
   #include <GL/glew.h>  <br/>
   to Viewer.h and Shader.h


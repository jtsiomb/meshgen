meshgen - a procedural mesh generation tool
===========================================

About
-----
Meshgen is a tool to simplify the creation of procedural meshes. It uses the
system C++ compiler, to build a user-supplied snippet of C++ code, which needs
to expose a `generate` function. The user-supplied function is called whenever
the user hits the `R` key, and the generated mesh is rendered in the viewport.
The current mesh is saved to a Wavefront OBJ file when pressing the `W` key.
Any generated textures are saved as individual image files, and referenced in
the exported OBJ material file.

![shot](http://nuclear.mutantstargoat.com/sw/meshgen/img/meshgen640.jpg)

Please note that this tool is meant as an aid for creating procedural geometry
with my `Mesh` class and mesh generation utility functions (see `mesh/`
directory in http://github.com/jtsiomb/dropcode). It might seem awkward to
require the use of a seemingly random collection of framework code in the
generator functions rather than a purpose-designed mesh generation API, but
that's by design for the intended purpose of this tool.

To learn how to use this tool, see the example mesh generator `test.cc` in the
root directory, and consult the header files `mesh.h` and `meshgen.h` for the
available mesh generation and manipulation functions.

All the vector/matrix math types and operations are provided by my current math
library: http://github.com/jtsiomb/gph-math. So you need to install that first.

License
-------
Copyright (C) 2018 John Tsiombikas <nuclear@member.fsf.org>

Meshgen is free software, you may use, copy, modify, and redistribute it under
the terms of the GNU General Public License v3, or any later version published
by the Free Software Foundation. See COPYING for detail.

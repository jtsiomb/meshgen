meshgen - a procedural mesh generation tool
===========================================

About
-----
Meshgen is a tool to simplify the creation of procedural meshes. It uses the
system C++ compiler, to build a user-supplied snippet of C++ code,
which needs to expose a `generate` function. The user-supplied function is
called whenever the user hits the `R` key, and the generated mesh is rendered in
the viewport. The current mesh is saved to a Wavefront OBJ file when pressing
the `W` key.

License
-------
Meshgen is free software, you may use, copy, modify, and redistribute it under
the terms of the GNU General Public License v3, or any later version published
by the Free Software Foundation. See COPYING for detail.

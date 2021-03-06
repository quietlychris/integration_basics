## Integration Basics in Rust

This repository contains code that has been ported to Rust from C++ based on Glenn Fiedler's Gaffer On Games "Integration Basics" article (see [here](https://gafferongames.com/post/integration_basics/)). I found it to be highly instructive in helping to transition from a basic understanding of Euler's method to a better understanding of RK4 and similar numerical methods. 

The initial compilation time may take a while, since the `inline-python` dependency that is being used for embedding Python's matplotlib has a fairly complicated compilation step. Subsequent iterations (i.e. each time you make an edit after the intial compilation) should be much faster. I'd recommend enabling the cargo `--release` flag on your first build, since it will make the for-loops run much more quickly. 

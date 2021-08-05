# MPI-OpenMP-Parallel-Mandelbrot-Set
Parallel Mandelbrot Set Algorithm implemented in C and C++ using MPI and OpenMP

A repeating calculation is performed for each pixel (each x and y points) in the area. Depending on the behaviour of that calculation, a color is chosen for that pixel. The algorithm itself is based on the number of iterations necessary to determine whether the orbit sequence tends to infinity or not. In my case I was using 255, and later 512 iterations to test my results.

Since a complicated calculation must be performed on every single pixel, total running time can easily become huge for a serial program. To increase efficiency of our program we are going to use MPI to be able to split tasks between nodes on the supercomputer and openMP to create threads and parallelize work on every of those nodes.

# Results:
![Mandelbrot Set Image](./Mandelbrot.jpg?raw=true "Mandelbrot Set")

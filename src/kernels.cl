__kernel void add(__global float* buffer, float scalar) {
    buffer[get_global_id(0)] += scalar;
}

__kernel void addtwo(__global float* buffer, float scalar) {
    buffer[get_global_id(0)] += scalar;
}

__kernel void gauss_seidel_step(__global float* buffer_new, __global float* buffer_old, float factor_a, float factor_c) {
    uint x = get_global_id(0);
    uint y = get_global_id(1);
    uint n_x = get_global_size(0);
    uint n_y = get_global_size(1);

    // Linear solving is impossible on buffer edges
    if (x != 0 && y != 0 && x != n_x && y != n_y) {
        buffer_new[x, y] = (buffer_old[x, y] + 
        factor_a * (buffer_new[x-1, y] + buffer_new[x+1, y] + buffer_new[x, y-1] + buffer_new[x, y+1])) 
        / factor_c;
    }
}
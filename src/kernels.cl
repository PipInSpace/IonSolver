__kernel void add(__global float* buffer, float scalar) {
    uint index = get_global_id(0) + (get_global_id(1) * get_global_size(1));
    buffer[index] += scalar;
}

__kernel void addtwo(__global float* buffer, float scalar) {
    uint index = get_global_id(0) + (get_global_id(1) * get_global_size(1));
    buffer[index] += scalar;
}

__kernel void gauss_seidel_step(__global float* buffer_new, __global float* buffer_old, float factor_a, float factor_c) {
    uint x = get_global_id(0);
    uint y = get_global_id(1);
    uint n_x = get_global_size(0);
    uint n_y = get_global_size(1);

    // Linear solving is impossible on buffer edges
    if (x != 0 && y != 0 && x != n_x-1 && y != n_y-1) {
        buffer_new[x + y*n_y] = (buffer_old[x + y*n_y] + 
        factor_a * (buffer_new[x - 1 + y*n_y] + buffer_new[x + 1 + y*n_y] + buffer_new[x + (y-1)*n_y] + buffer_new[x + (y+1)*n_y])) 
        / factor_c;
    }
}

__kernel void advect(int b, __global float* buffer_a, __global float* buffer_b, __global float* force_x, __global float* force_y, float dt) {
    uint x_index = get_global_id(0);
    uint y_index = get_global_id(1);
    uint n_x = get_global_size(0);
    uint n_y = get_global_size(1);

    // Advection is impossible at buffer edges
    // Valid x/y value indexes range from 1 to (n_x/n_y-2)
    if (x_index != 0 && y_index != 0 && x_index != n_x-1 && y_index != n_y-1) {
        uint buffer_index = x_index + (y_index * n_y);

        float x = as_float(x_index) - dt * force_x[buffer_index];
        float y = as_float(y_index) - dt * force_y[buffer_index];

        if (x < 0.5f) {
            x = 0.5f;
        }
        if (x > as_float(n_x-2) + 0.5f) {
            x = as_float(n_x-2) + 0.5f;
        }
        uint x_idx0 = as_uint(x);
        uint x_idx1 = x_idx0 + 1;

        if (y < 0.5f) {
            y = 0.5f;
        }
        if (y > as_float(n_y-2) + 0.5f) {
            y = as_float(n_y-2) + 0.5f;
        }
        uint y_idx0 = as_uint(y);
        uint y_idx1 = y_idx0 + 1;

        float fac_s1 = x - as_float(x_idx0);
        float fac_s0 = 1.0f - fac_s1;
        float fac_t1 = y - as_float(y_idx0);
        float fac_t0 = 1.0f - fac_t1;

        buffer_a[buffer_index] = fac_s0 * (fac_t0 * buffer_b[x_idx0 + y_idx0*n_y] + fac_t1 * buffer_b[x_idx0 + y_idx1*n_y])
        + fac_s1 * (fac_t0 * buffer_b[x_idx1 + y_idx0*n_y] + fac_t1 * buffer_b[x_idx1 + y_idx1*n_y]);
    }
}
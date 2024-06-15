#include "dgemm.h"
#include <cstdio>
#include <cstdlib>
#include <immintrin.h>

void dgemm(float alpha, const float *a, const float *b, float beta, float *c) {
    __m256 alpha_vec = _mm256_set1_ps(alpha);
    
    
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            float partial_sum_array[8] = {0, 0, 0, 0, 0,0,0,0};
            __m256 partial_sum = _mm256_set1_ps(0);
            c[i * MATRIX_SIZE + j] *= beta;
            for (int k = 0; k < MATRIX_SIZE - 8; k+=8) {
                __m256 a_vec = _mm256_loadu_ps(a+(i*MATRIX_SIZE)+k);
                __m256 b_vec = _mm256_loadu_ps(b+(j*MATRIX_SIZE)+k);
                __m256 mul = _mm256_mul_ps(a_vec,b_vec);
                __m256 mul1 = _mm256_mul_ps(alpha_vec,mul);
                partial_sum = _mm256_add_ps(mul1,partial_sum);
                for(int l = 0; l < 8; l++) {
                    c[i * MATRIX_SIZE + j] += alpha * a[i * MATRIX_SIZE + k+l] * b[j * MATRIX_SIZE + k+l];
                }
                //c[i * MATRIX_SIZE + j] += alpha * a[i * MATRIX_SIZE + k] * b[j * MATRIX_SIZE + k];
            }
            _mm256_store_ps(partial_sum_array,partial_sum);
            for (int l = 0; l < 8; l++) {
                c[i * MATRIX_SIZE + j] += partial_sum_array[l];
            }
            for (int k = MATRIX_SIZE - (MATRIX_SIZE % 8); k < MATRIX_SIZE; k++) {
                c[i * MATRIX_SIZE + j] += alpha * a[i * MATRIX_SIZE + k] * b[j * MATRIX_SIZE + k];
            }
        }
    }
}

int main(int, char **) {
    float alpha, beta;

    // mem allocations
    int mem_size = MATRIX_SIZE * MATRIX_SIZE * sizeof(float);
    auto a = (float *) aligned_alloc(16,mem_size);
    auto b = (float *) aligned_alloc(16,mem_size);
    auto c = (float *) aligned_alloc(16,mem_size);

    // check if allocated
    if (nullptr == a || nullptr == b || nullptr == c) {
        printf("Memory allocation failed\n");
        if (nullptr != a) free(a);
        if (nullptr != b) free(b);
        if (nullptr != c) free(c);
        return 0;
    }

    generateProblemFromInput(alpha, a, b, beta, c);

    std::cerr << "Launching dgemm step." << std::endl;
    // matrix-multiplication
    dgemm(alpha, a, b, beta, c);

    outputSolution(c);

    free(a);
    free(b);
    free(c);
    return 0;
}

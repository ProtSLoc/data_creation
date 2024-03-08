#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <cuda_runtime.h>

// CUDA kernel to calculate sequence identity on GPU
__global__ void calculate_sequence_identity_gpu(const char* seq1, const char* seq2, float* result, size_t seq_len) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < seq_len) {
        result[idx] = (seq1[idx] == seq2[idx]) ? 1.0f : 0.0f;
    }
}

// Function to calculate sequence identity using GPU
float calculate_sequence_identity_cuda(const char* seq1, const char* seq2, size_t seq_len) {
    const int threads_per_block = 256;
    const int blocks = (seq_len + threads_per_block - 1) / threads_per_block;

    float* d_result;
    cudaMalloc((void**)&d_result, seq_len * sizeof(float));
    calculate_sequence_identity_gpu<<<blocks, threads_per_block>>>(seq1, seq2, d_result, seq_len);

    float* h_result = new float[seq_len];
    cudaMemcpy(h_result, d_result, seq_len * sizeof(float), cudaMemcpyDeviceToHost);

    float sum = 0.0f;
    for (size_t i = 0; i < seq_len; ++i) {
        sum += h_result[i];
    }

    delete[] h_result;
    cudaFree(d_result);

    return (100.0f * sum) / (seq_len / 2.0f);
}

// Function to check sequence validity
bool check_sequence_validity_cuda(const char* seq1, const std::vector<std::string>& dataset2_sequences, float threshold) {
    for (const auto& seq2 : dataset2_sequences) {
        float score = calculate_sequence_identity_cuda(seq1, seq2.c_str(), std::min(seq1.length(), seq2.length()));
        if (score > threshold) {
            return false;
        }
    }
    return true;
}

// Function to process a sequence
bool process_sequence_cuda(int index, const std::string& row, const std::vector<std::string>& dataset2_sequences, float threshold) {
    return check_sequence_validity_cuda(row.c_str(), dataset2_sequences, threshold);
}

int main() {
    // Load protein sequences from CSV files
    std::vector<std::string> dataset1_sequences = load_sequences_from_csv("/content/data_creation/Dataset/deeploc_data.csv");
    std::vector<std::string> dataset2_sequences = load_sequences_from_csv("/content/data_creation/Dataset/hpa_testset.csv");

    float threshold = 30.0;  // Set your desired threshold

    std::vector<bool> results(dataset1_sequences.size(), false);

    // Process sequences in parallel using threads
    #pragma omp parallel for
    for (size_t i = 0; i < dataset1_sequences.size(); ++i) {
        results[i] = process_sequence_cuda(i, dataset1_sequences[i], dataset2_sequences, threshold);
    }

    // Print results
    for (size_t i = 0; i < dataset1_sequences.size(); ++i) {
        std::cout << "Sequence " << i << ": " << (results[i] ? "Valid" : "Not Valid") << std::endl;
    }

    return 0;
}

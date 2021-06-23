#ifndef POINT_PROCESS_H
#define POINT_PROCESS_H

#include <algorithm>
#include <vector>
#include <random>
#include <string>

enum PointProcessType{RSA, Random, File};

void binomial_point_process(std::vector<Vec2d>& points, unsigned int N, unsigned int seed) {
    std::mt19937 generator;
    generator.seed(seed);
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    points.clear();
    for (unsigned int s = 0; s < N; s++) {
        points.push_back(Vec2d(distribution(generator), distribution(generator)));
    }
}

void poisson_point_process(std::vector<Vec2d>& points, double intensity, double area, unsigned int seed) {
    std::mt19937 generator;
    generator.seed(seed);
    double mean = intensity * area;
    std::poisson_distribution<int> distribution(mean);
    int N = distribution(generator);
    binomial_point_process(points, N, seed);
}

void rsa_point_process(std::vector<Vec2d>& points, double max_distance, bool is_periodic, unsigned int seed, unsigned int max_trials) {
    points.clear();
    const int cell_id[3] = {0, -1, 1};
    int size_cell_id = 1;
    if (is_periodic) size_cell_id = 3;
    std::mt19937 generator;
    generator.seed(seed);
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    const double max_distance_squared = max_distance * max_distance;
    for (unsigned int i = 0; i < max_trials; i++) {
        Vec2d candidate_point = Vec2d(distribution(generator), distribution(generator));
        bool reject = false;
        for (unsigned int s = 0; s < points.size() && !reject; s++) {
            for (int j = 0; j < size_cell_id && !reject; j++) {
                int rx = cell_id[j];
                for (int k = 0; k < size_cell_id && !reject; k++) {
                    int ry = cell_id[k];
                    Vec2d point = points[s] + Vec2d(rx, ry);
                    double distance_squared = (point - candidate_point).length_squared();
                    if (distance_squared < max_distance_squared) {
                        reject = true;
                    }
                }
            }
        }

        if (!reject) {
            points.push_back(candidate_point);
        }
    }
}

void file_point_process(std::vector<Vec2d>& points, const std::string& filename) {
    std::ifstream file(filename.c_str(), std::ios::binary);
    float coord;
    unsigned int i = 0;
    Vec2d point(0, 0);
    while (file.read(reinterpret_cast<char*>(&coord), sizeof(float))) {
        point[i] = coord;
        if (i == 1) {
            points.push_back(point);
        }
        i = (i + 1) % 2;
    }
}

#endif

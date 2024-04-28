#ifndef GENETICALGORITHM_DATAPARSERQUICK_H
#define GENETICALGORITHM_DATAPARSERQUICK_H
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include "Individuals.h"

std::vector<bool> parseBoolVector(const std::string& line, size_t n) {
    std::vector<bool> result;
    result.reserve(n);
    for (size_t i = 0; i < line.size(); i += 2) {
        result.push_back(line[i] == '1');
    }
    return result;
}

std::vector<double> parseDoubleVector(const std::string& line, size_t n) {
    std::vector<double> result;
    result.reserve(n);
    std::istringstream iss(line);
    double value;
    char comma;
    while (iss >> value) {
        result.push_back(value);
        iss >> comma;
    }
    return result;
}

void formAdjacencyMatrix(std::vector<std::vector<bool>> &matrix) {
    for (int i = 0; i < matrix.size(); ++i)
        for (int j = i; j < matrix.size(); ++j) {
            matrix[j][i] = 1 - matrix[i][j];
            if (i == j)
                continue;
            matrix[i][j] = 1 - matrix[i][j];
        }
}

std::pair<std::vector<std::vector<bool>>, std::vector<double>> parseCSV(const std::string& dataPath) {
    std::ifstream file(dataPath);
    if (!file.is_open()) {
        std::cerr << "Error opening file" << std::endl;
        return {};
    }

    std::vector<std::vector<bool>> boolMatrix;
    std::vector<double> lastLineDoubles;
    std::string line;

    if (std::getline(file, line)) {
        size_t n = std::count(line.begin(), line.end(), ',') + 1;
        boolMatrix.push_back(parseBoolVector(line, n));

        for (auto _ = 0; _ < n - 1; ++_) {
            std::getline(file, line);
            if (file.peek() != EOF) {
                boolMatrix.push_back(parseBoolVector(line, n));
            }
        }
        std::getline(file, line);
        lastLineDoubles = parseDoubleVector(line, n); // Parse the last line
    }

    return {boolMatrix, lastLineDoubles};
}

void writeCSV(/*const std::vector<std::vector<bool>>& matrix, */
              const std::vector<Individual>& matrix, const std::string& dataPath, int n) {

    FILE* file = fopen(dataPath.c_str(), "w");

    if (file == nullptr) {
        std::cerr << "Error opening file for writing." << std::endl;
        return;
    }

    for (auto i = 1; i < n + 1; ++i){
        fprintf(file, ",TH%d", i);
    }
    fprintf(file, ",sum(w)\n");

    for (int i = matrix.size() - 1; i >= 0; --i) {

        fprintf(file, "GH%d", matrix.size() - i);
        auto end_itr = matrix[i].nodes.end();
        auto itr = matrix[i].nodes.begin();

        for (auto c = 0; c < n; ++c){
            if (itr != end_itr and c == *itr) {
                fprintf(file, ",1");
                ++itr;
            }
            else
                fprintf(file, ",0");
        }
        fprintf(file, ",%.3f\n", matrix[i].getWeight());
    }

    fclose(file);
}

#endif //GENETICALGORITHM_DATAPARSERQUICK_H

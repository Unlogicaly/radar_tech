#include <iostream>
#include <chrono>
#include "DataParserQuick.h"
#include "GeneticAlgorithm.h"

int main(int argc, char *argv[]) {

    auto start = std::chrono::high_resolution_clock::now();

    const char* dataPath = argv[1];

    std::cout << "datapath " << dataPath << std::endl;//"C:/LOGS/input_with_weights.csv"; // replace this with the actual path
    auto [matrix, weights] = parseCSV(dataPath);
    formAdjacencyMatrix(matrix);

    std::cout << "input time, ms: " << std::chrono::duration_cast<std::chrono::microseconds>
        (std::chrono::high_resolution_clock::now() - start).count()/1000 << std::endl;

    start = std::chrono::high_resolution_clock::now();

    Genetics genetics{};
    auto result = genetics.geneticAlgorithm(matrix, weights);

    std::cout << "algorithm time, ms: " << std::chrono::duration_cast<std::chrono::microseconds>
            (std::chrono::high_resolution_clock::now() - start).count()/1000 << std::endl;

    start = std::chrono::high_resolution_clock::now();

    std::sort(result.begin(), result.end(), [](const Individual &lhs, const Individual &rhs) {
        return lhs.getWeight() < rhs.getWeight();
    });

    writeCSV(result, argv[2], weights.size());
    for (auto &ind: result) {
        std::cout << "nodes: {";
        for (auto &node: ind.nodes)
            std::cout << node << ", ";
        std::cout << "}, weight: " << ind.getWeight() << "\n";
    }

    std::cout << "output time, ms: " << std::chrono::duration_cast<std::chrono::microseconds>
            (std::chrono::high_resolution_clock::now() - start).count()/1000 << std::endl;

    return 0;
}

#include "DataParserQuick.h"
#include "GeneticAlgorithm.h"

int main(int argc, char *argv[]) {

    const char* dataPath = argv[1];
    auto [matrix, weights] = parseCSV(dataPath);
    formAdjacencyMatrix(matrix);

    Genetics genetics{};
    auto result = genetics.geneticAlgorithm(matrix, weights);

    std::sort(result.begin(), result.end(), [](const Individual &lhs, const Individual &rhs) {
        return lhs.getWeight() < rhs.getWeight();
    });

    writeCSV(result, argv[2], weights.size());

    return 0;
}

#ifndef GENETICALGORITHM_POPULATION_H
#define GENETICALGORITHM_POPULATION_H

#include "Individuals.h"
#include <vector>

class Population {
public:
    Population(const std::vector<Individual> &individuals);

    void addIndividual(const Individual& individual);

    int size() const { return individuals.size(); }

    const Individual& operator [] (int i) const { return individuals[i]; }

    std::vector<double> getWeights() const { return weights; };

private:
    std::vector<Individual> individuals;
    std::vector<double> weights;
    double total_weight;
};

#endif //GENETICALGORITHM_POPULATION_H

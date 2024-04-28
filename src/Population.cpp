#include "Population.h"

Population::Population(const std::vector <Individual> &individuals)  {
    this->individuals = individuals;

    int n_individuals = individuals.size();
    weights = std::vector<double>(n_individuals);
    total_weight = 0;

    for (int i = 0; i < n_individuals; ++i) {
        weights[i] = individuals[i].getWeight();
        total_weight += weights[i];
    }
}

void Population::addIndividual(const Individual &individual) {
    individuals.push_back(individual);
    total_weight += individual.getWeight();
    weights.push_back(individual.getWeight());
}
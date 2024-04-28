#ifndef GENETICALGORITHM_GENETICALGORITHM_H
#define GENETICALGORITHM_GENETICALGORITHM_H

#include <vector>
#include <set>
#include <algorithm>
#include <random>
#include "ctime"
#include "Individuals.h"
#include "Population.h"


class Genetics {
public:
    Genetics() {
        normalDistribution = std::normal_distribution<double>{0, 1};
        g = std::mt19937 (rd());
    }

    std::vector<Individual> geneticAlgorithm(std::vector<std::vector<bool>> &matrix,
                                             std::vector<double> &weights);

private:
    std::random_device rd;
    std::mt19937 g;
    std::normal_distribution<double> normalDistribution;
    std::set<double> top;

    Population selectionByWeight(const Population& population);

    std::vector<std::pair<Individual, Individual>> parentSelectionOut(const Population& population, double skip_parent_prob=0.);

    int matrix_sum(const std::vector<std::vector<bool>>& matrix) {
        int sum = 0;
        for (int i = 0; i < matrix.size(); ++i) {
            sum += std::accumulate(matrix[i].begin(), matrix[i].end(), 0);
        }
        return sum;
    }

    int argmax(const std::vector<double>& v) {
        auto max = max_element(v.begin(), v.end());
        int argmaxVal = distance(v.begin(), max);
        return argmaxVal;
    }

    std::vector<int> arange(int n) {
        std::vector<int> out(n);
        for (int i = 0; i < n; ++i)
            out[i] = i;
        return out;
    }

    void erase_column(std::vector<std::vector<bool>>& matrix, int n) {
        for (int i = 0; i < matrix.size(); ++i)
            matrix[i].erase(matrix[i].begin() + n);
    }
    void erase_row(std::vector<std::vector<bool>>& matrix, int n) {
        matrix.erase(matrix.begin() + n);
    }

    void updateSet(std::set<int> & to_update,
                   const std::vector<bool>& vec) {
        std::vector<int> idx_false{};
        idx_false.reserve(vec.size());

        for (int i = 0; i < vec.size(); ++i) {
            if (vec[i])
                idx_false.push_back(i);
        }

        to_update.insert(idx_false.begin(), idx_false.end());
    }

    std::vector<int> randomChoiceReplaceFalse(const std::vector<int>& nodes, int sizeOfShuffle) {
        auto idx = arange(sizeOfShuffle);
        std::shuffle(idx.begin(), idx.end(), g);

        std::vector<int> out(idx.size());
        for (int i = 0; i < idx.size(); ++i) {
            out[i] = nodes[idx[i]];
        }
        return out;
    }

    void updateCandidates(std::vector<int>& candidates,
                                   const std::vector<std::vector<bool>>& matrix,
                                   const std::vector<double>& weights,
                                   std::vector<int>& indexes);

    std::vector<int> truncationGreedy(std::vector<std::vector<bool>> matrix,
                                      std::vector<double>& weights);


    std::vector<int> expansionGreedy(std::vector<int> &indexes, const std::vector<double> &weights,
                                     const std::vector<std::vector<bool>> &matrix);

    std::vector<int> mutation(std::vector<int>& nodes_list,
                              int max_node_idx, double prob=0.1,
                              double strength=0.1);

    Individual breed(const Individual& lhs, const Individual& rhs,
                     const std::vector<std::vector<bool>>& matrix,
                     const std::vector<double>& weights);


    Individual createIndividual(std::vector<std::vector<bool>> matrix,
                                std::vector<double>& weights);

    Population createPopulation(std::vector<std::vector<bool>>& matrix,
                                std::vector<double>& weights,
                                int population_size=20);

    std::vector<Individual> breedAllParents(const std::vector<std::pair<Individual, Individual>>& parents,
                                            const std::vector<std::vector<bool>>& matrix,
                                            const std::vector<double>& weights);

    std::vector<int> truncationRandom(std::vector<std::vector<bool>> matrix,
                                      std::vector<double>& weights,
                                      double temperature=2.25, double alpha=0.5);

    std::vector<int> expansionRandom(std::vector<int> &indexes,
                                      const std::vector<double> &weights,
                                      const std::vector<std::vector<bool>> &matrix,
                                      double temperature=2.25, double alpha=0.5);



};
#endif //GENETICALGORITHM_GENETICALGORITHM_H

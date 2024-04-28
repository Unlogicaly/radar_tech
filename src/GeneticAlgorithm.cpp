#include <iostream>
#include <chrono>
#include "GeneticAlgorithm.h"
#include <omp.h>

Population Genetics::selectionByWeight(const Population &population)  {
    auto weights = population.getWeights();

    std::vector<int> indexes(population.size());
    for (int i = 0; i < indexes.size(); ++i)
        indexes[i] = i;

    std::shuffle(indexes.begin(), indexes.end(), g);

    std::vector<Individual> individuals(indexes.size()/2);
    for (int i = 0; i < indexes.size()/2; ++i) {
        if (weights[i] > weights[weights.size() - i - 1])
            individuals[i] = population[i];
        else
            individuals[i] = population[weights.size() - i - 1];
    }

    auto result = Population{individuals};

    return result;
}

std::vector<std::pair<Individual, Individual>>
Genetics::parentSelectionOut(const Population &population, double skip_parent_prob)  {
    std::uniform_int_distribution<int> dist(0, population.size() - 1);



    std::vector<int> first_parent_idx(population.size()/2);

    for (int i = 0; i < population.size()/2; ++i)
        first_parent_idx[i] = dist(g);

    std::vector<std::pair<Individual, Individual>> result;
    result.reserve(first_parent_idx.size());

    clock_t start = clock();

    for (auto& i : first_parent_idx) {
        int max_diff = 0;
        int max_diff_idx = 0;

        for (int second_parent_idx = 0; second_parent_idx < population.size(); ++second_parent_idx) {
            if (i == second_parent_idx)
                continue;

            if (normalDistribution(g) < skip_parent_prob)
                continue;

            std::vector<int> diff;
            std::set_symmetric_difference(population[i].nodes.begin(),
                                          population[i].nodes.end(),
                                          population[second_parent_idx].nodes.begin(),
                                          population[second_parent_idx].nodes.end(),
                                          std::back_inserter(diff));

            int len_diff = diff.size();
            if (len_diff > max_diff) {
                max_diff = len_diff;
                max_diff_idx = second_parent_idx;
            }

        }

        result.emplace_back(population[i], population[max_diff_idx]);
    }

    return result;
}

void Genetics::updateCandidates(std::vector<int> &candidates, const std::vector<std::vector<bool>> &matrix,
                                const std::vector<double> &weights, std::vector<int> &indexes)  {

    std::vector<double> column_sums(candidates.size());

    for (size_t i = 0; i < candidates.size(); ++i) {
        for (size_t j = 0; j < candidates.size(); ++j) {
            if (matrix[candidates[i]][candidates[j]]) {
                column_sums[j]++;
            }
        }
        column_sums[i] /= (weights[candidates[i]] + 1e-5);
    }

    int current_node_index = argmax(column_sums);
    int current_node = candidates[current_node_index];

    indexes.push_back(current_node);

    candidates.erase(
            std::remove_if(
                    candidates.begin(),
                    candidates.end(),
                    [&](const int& candidate) {
                        return matrix[current_node][candidate]; }
            ),
            candidates.end()
    );
}

std::vector<int> Genetics::mutation(std::vector<int> &nodes_list, int max_node_idx, double prob, double strength)  {

    if (normalDistribution(g) < prob)
        return nodes_list;

    std::uniform_int_distribution<int> dist(0, ceil(max_node_idx*strength - 1));
    int nodes_number_to_add = dist(g);
    int nodes_number_to_remove = dist(g);

    auto idx = arange(max_node_idx);

    std::vector<int> nodes_not_included{};
    std::set_difference(idx.begin(), idx.end(),
                        nodes_list.begin(), nodes_list.end(),
                        std::back_inserter(nodes_not_included));

    int n_nodes_to_add = std::min({nodes_number_to_add, (int)nodes_not_included.size()});
    auto nodes_to_add = randomChoiceReplaceFalse(nodes_not_included, n_nodes_to_add);

    nodes_list.insert(nodes_list.begin(), nodes_to_add.begin(), nodes_to_add.end());

    std::shuffle(nodes_list.begin(), nodes_list.end(), g);

    nodes_list.erase(nodes_list.begin() + nodes_number_to_remove, nodes_list.end());
    return nodes_list;
}

Individual Genetics::breed(const Individual &lhs, const Individual &rhs,
                           const std::vector<std::vector<bool>> &matrix,
                           const std::vector<double> &weights)  {

    auto start = std::chrono::high_resolution_clock::now();

    std::vector<int> all_nodes;
    std::set_union(lhs.nodes.begin(), lhs.nodes.end(),
                   rhs.nodes.begin(), rhs.nodes.end(),
                   std::back_inserter(all_nodes));

    all_nodes = mutation(all_nodes, matrix.size());

    int sub_size = all_nodes.size();
    std::vector<std::vector<bool>> sub_matrix(sub_size, std::vector<bool>(sub_size));
    std::vector<double> sub_weights(sub_size);
    for (int i = 0; i < sub_size; ++i) {
        for (int j = 0; j < sub_size; ++j)
            sub_matrix[i][j] = matrix[all_nodes[i]][all_nodes[i]];
        sub_weights[i] = weights[all_nodes[i]];
    }

    auto idxs = truncationRandom(sub_matrix, sub_weights);

    std::vector<int> truncated(idxs.size());
    for (int i = 0; i < idxs.size(); ++i)
        truncated[i] = all_nodes[idxs[i]];

    auto expanded_nodes = expansionRandom(truncated, weights, matrix);

    double weight_sum = 0;
    for (int expanded_node : expanded_nodes)
        weight_sum += weights[expanded_node];

    return Individual{std::set(expanded_nodes.begin(), expanded_nodes.end()), weight_sum};
}


std::vector<Individual> Genetics::breedAllParents(const std::vector<std::pair<Individual, Individual>>& parents,
                                        const std::vector<std::vector<bool>>& matrix,
                                        const std::vector<double>& weights) {
    std::vector<Individual> out(parents.size());

    #pragma omp parallel for
        for (auto i = 0; i < parents.size(); ++i) {
            out[i] = breed(parents[i].first, parents[i].second, matrix, weights);
        }
//    for (const auto & parent : parents) {
//        out.push_back(breed(parent.first, parent.second, matrix, weights));
//    }

    return out;
}

Individual Genetics::createIndividual(std::vector<std::vector<bool>> matrix, std::vector<double> &weights)  {
    std::vector<int> res;
    res.reserve(matrix.size());

    auto indexes = arange(matrix.size());

    while (indexes.size() > 0) {
        std::uniform_int_distribution<int> dist(0, indexes.size() - 1);
        int idx = dist(g);
        res.push_back(indexes[idx]);

        std::vector<int> where_matrix_1;
        where_matrix_1.reserve(matrix.size());

        for (int i = 0; i < matrix.size(); ++i) {
            if (matrix[idx][i] == 1)
                where_matrix_1.push_back(i);
        }

        for (int i = where_matrix_1.size() - 1; i >= 0; --i) {
            int j = where_matrix_1[i];
            indexes.erase(indexes.begin() + j);
            erase_column(matrix, j);
            erase_row(matrix, j);
        }
    }

    double sum = 0;
    for (int re : res) {
        sum += weights[re];
    }
    return Individual{std::set(res.begin(), res.end()), sum};
}

Population
Genetics::createPopulation(std::vector<std::vector<bool>> &matrix, std::vector<double> &weights, int population_size) {

    std::vector<Individual> inds{};
    for (int i = 0; i < population_size; ++i)
        inds.push_back(createIndividual(matrix, weights));

    return Population{inds};
}



class IndCompHeap{
public:
    bool operator () (const Individual &lhs, const Individual &rhs) {
        return lhs.getWeight() > rhs.getWeight();
    }
};

std::vector<Individual> Genetics::geneticAlgorithm(std::vector<std::vector<bool>> &matrix, std::vector<double> &weights) {
    auto population = createPopulation(matrix, weights, 50);

    int no_change = 0;
    // todo config
    int no_change_config = 15;
    int config_keep_best = 5;

    std::set<std::set<int>> unique_bodies;
    std::vector<Individual> top_scores_heap{};


    while (no_change < no_change_config) {

        auto selected_population = selectionByWeight(population);

        auto parents = parentSelectionOut(population);

        auto children = breedAllParents(parents, matrix, weights);

        for (auto& child: children)
            selected_population.addIndividual(child);

        population = selected_population;

        bool changed = false;

        for(int i = 0; i < population.size(); ++i) {

            if (unique_bodies.find(population[i].nodes) == unique_bodies.end()) {
                unique_bodies.insert(population[i].nodes);
                top_scores_heap.push_back(population[i]);
                std::push_heap(top_scores_heap.begin(), top_scores_heap.end(), IndCompHeap());
                if (top_scores_heap.size() > config_keep_best) {
                    std::pop_heap(top_scores_heap.begin(), top_scores_heap.end(), IndCompHeap());
                    if (!(top_scores_heap[top_scores_heap.size() - 1] == population[i]))
                        changed = true;
                    top_scores_heap.pop_back();
                }
            }
        }

        if (changed)
            no_change = 0;
        else
            no_change += 1;
    }
    return top_scores_heap;
}

std::vector<int>
Genetics::truncationRandom(std::vector<std::vector<bool>> matrix, std::vector<double> &weights, double temperature,
                           double alpha) {
    auto matrix_size = matrix.size();

    for (int i = 0; i < matrix_size; ++i) {
        matrix[i][i] = false;
    }

    std::vector<int> deleted_indexes;
    deleted_indexes.reserve(matrix_size);

    auto indexes = arange(matrix_size);

    std::vector<double> coefs(matrix_size);

    while (matrix_sum(matrix) > 0) {

        for (int i = 0; i < matrix_size; ++i) {
            coefs[i] = std::accumulate(matrix[i].begin(), matrix[i].end(), 0);
            if (coefs[i] < 1e-2)
                coefs[i] = -1e6;
            else
                coefs[i] -= alpha*weights[i];

            coefs[i] = std::exp(coefs[i] / temperature);
        }

        double coef_sum = std::accumulate(coefs.begin(), coefs.end(), 0.);
        for (auto i = 0; i < matrix_size; ++i) {
            coefs[i] /= coef_sum;
        }

        std::discrete_distribution<> d(coefs.begin(), coefs.end());
        auto current_node = indexes[d(g)];
        for (int i = 0; i < matrix.size(); ++i) {
            matrix[current_node][i] = false;
            matrix[i][current_node] = false;
        }
        deleted_indexes.push_back(current_node);
    }

    std::vector<int> idx(matrix.size());
    for (int i = 0; i < matrix.size(); ++i)
        idx[i] = i;

    std::vector<int> res;

    std::sort(deleted_indexes.begin(), deleted_indexes.end());
    std::set_difference(idx.begin(), idx.end(),
                        deleted_indexes.begin(), deleted_indexes.end(),
                        std::back_inserter(res));
    return res;
}


std::vector<int> Genetics::truncationGreedy(std::vector<std::vector<bool>> matrix, std::vector<double> &weights) {

    for (int i = 0; i < matrix.size(); ++i) {
        matrix[i][i] = false;
    }

    std::vector<int> deleted_indexes;
    deleted_indexes.reserve(matrix.size());

    while (matrix_sum(matrix) > 0) {

        std::vector<double> sum_vec(matrix.size());

        for (int i = 0; i < matrix.size(); ++i) {
            sum_vec[i] = std::accumulate(matrix[i].begin(), matrix[i].end(), 0);
            sum_vec[i] /= (weights[i] + 1e-5);
        }

        int current_node = argmax(sum_vec);
        for (int i = 0; i < matrix.size(); ++i) {
            matrix[current_node][i] = false;
            matrix[i][current_node] = false;
        }
        deleted_indexes.push_back(current_node);
    }

    std::vector<int> idx(matrix.size());
    for (int i = 0; i < matrix.size(); ++i)
        idx[i] = i;

    std::vector<int> res;

    std::sort(deleted_indexes.begin(), deleted_indexes.end());
    std::set_difference(idx.begin(), idx.end(),
                        deleted_indexes.begin(), deleted_indexes.end(),
                        std::back_inserter(res));
    return res;
}

std::vector<int> Genetics::expansionRandom(std::vector<int> &indexes, const std::vector<double> &weights,
                                           const std::vector<std::vector<bool>> &matrix, double temperature,
                                           double alpha) {
    auto all_indexes = arange(matrix.size());
    std::set<int> indexes_to_delete{};

    for (auto& idx: indexes) {
        updateSet(indexes_to_delete, matrix[idx]);
    }

    std::vector<int> candidates{};

    std::set_difference(all_indexes.begin(), all_indexes.end(),
                        indexes_to_delete.begin(), indexes_to_delete.end(),
                        std::back_inserter(candidates));

    while (candidates.size() > 0) {
        std::vector<double> coefs(candidates.size());
        for (auto i = 0; i < candidates.size(); ++i) {
            for (auto j = 0; j < candidates.size(); ++j){
                coefs[i] -= double(matrix[candidates[i]][candidates[j]]);
            }

            coefs[i] += alpha*weights[candidates[i]];
            coefs[i] = std::exp(coefs[i] / temperature);
        }

        auto coef_sum = std::accumulate(coefs.begin(), coefs.end(), 0.);
        for (auto i = 0; i < candidates.size(); ++i) {
            coefs[i] /= coef_sum;
        }

        std::discrete_distribution<> d(coefs.begin(), coefs.end());
        auto current_node = candidates[d(g)];

        indexes.push_back(current_node);

        candidates.erase(
                std::remove_if(
                        candidates.begin(),
                        candidates.end(),
                        [&](const int& candidate) {
                            return matrix[current_node][candidate]; }
                ),
                candidates.end()
        );
    }

    return indexes;
}

std::vector<int>
Genetics::expansionGreedy(std::vector<int> &indexes, const std::vector<double> &weights,
                          const std::vector<std::vector<bool>> &matrix) {

    auto candidates = arange(matrix.size());
    std::set<int> indexes_to_delete{};

    for (auto& idx: indexes) {
        updateSet(indexes_to_delete, matrix[idx]);
    }

    std::set_difference(candidates.begin(), candidates.end(),
                        indexes_to_delete.begin(), indexes_to_delete.end(),
                        candidates.begin());

    while (candidates.size() > 0) {
        updateCandidates(candidates, matrix, weights, indexes);
    }

    return indexes;
}



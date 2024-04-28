#ifndef GENETICALGORITHM_INDIVIDUALS_H
#define GENETICALGORITHM_INDIVIDUALS_H

#include <set>

class Individual {
public:

    Individual() {
        this->nodes = {};
        this->weight = 0;
    }

    Individual(const std::set<int>& nodes,
               double weight) {
        this->nodes = nodes;
        this->weight = weight;
    }

    bool operator < (const Individual& other) const {
        return this->weight < other.weight;
    }

    bool operator == (const Individual& other) const {
        return this->nodes == other.nodes;
    }

    Individual& operator = (const Individual& other) {
        if (this == &other)
            return *this;
        this->nodes = other.nodes;
        this->weight = other.weight;

        return *this;
    }

    std::set<int> nodes;

    double getWeight() const {return weight; };
private:
    double weight;
};


#endif //GENETICALGORITHM_INDIVIDUALS_H

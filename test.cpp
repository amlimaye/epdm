#include "types.hpp"
#include "constants.hpp"
#include "utilities.hpp"
#include <iostream>
#include <numeric>
#include <random>
#include <cmath>

void add_population_to_ensemble(ensemble_t& ensemble, const population_t& population) {
    //make tuples representing all _new_ bimolecular reactions with this species
    std::list<std::tuple<population_t, population_t>> all_pairs;
    for (auto& existing_population : ensemble.populations) {
        all_pairs.push_back(std::make_tuple(existing_population,population));
    }
    all_pairs.push_back(std::make_tuple(population,population));

    //create relations for each pair
    for (auto& pair : all_pairs) {
        //unpack the tuple
        population_t p1,p2;
        std::tie(p1,p2) = pair;

        //make a relation_t object held by the _first_ member of the pair
        relation_t new_relation;
        new_relation.owner_population_ptr = &p1;

        //add reactions for this species pair
        new_relation.reactions = rxn_utilities::get_reactions_for_species_pair(p1,p2);
        new_relation.tot_partial_propensity = 0.0;
        for (auto& rxn: new_relation.reactions) {
            new_relation.tot_partial_propensity += rxn.partial_propensity;
        }

        //make a relation_address_t held by the _second_ member of the pair
        relation_address_t new_relation_address;
        new_relation_address.owner_population_ptr = &p2;
        new_relation_address.relation_ptr = &new_relation;

        //put them in the right places
        p1.relations.push_back(new_relation);
        p2.relation_addresses.push_back(new_relation_address);

        //update the propensities of p1
        p1.tot_partial_propensities.push_back(new_relation.tot_partial_propensity);
        p1.tot_propensity = std::accumulate(p1.tot_partial_propensities.begin(),
                                            p1.tot_partial_propensities.end(),
                                            0.0);

        //update the propensities of the ensemble
        ensemble.total_propensity += p1.tot_propensity;
    }

    //finally, add the population to the ensemble
    ensemble.populations.push_back(population);
}

void delete_population_from_ensemble(ensemble_t& ensemble) {
}

void take_timestep(ensemble_t& ensemble) {
    //sample the next time
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    double rand = distribution(ensemble.generator);
    double next_time = (1.0/ensemble.total_propensity)*std::log((1.0/rand));
}

ensemble_t initialize_ensemble() {
    //start at time zero with zero propensity
    ensemble_t ensemble;
    ensemble.current_time = 0.0;
    ensemble.total_propensity = 0.0;

    //ensembles always start with the void population
    add_population_to_ensemble(ensemble,pop_utilities::make_void_population());

    for (auto& elem : ensemble.populations) {
        pop_utilities::print_population(elem);
    }

    return ensemble;
}

int main() {
    auto ensemble = initialize_ensemble();
}

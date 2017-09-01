#include "types.hpp"
#include "constants.hpp"
#include "utilities.hpp"
#include <iostream>

std::list<reaction_t> get_reactions_for_species_pair(const population_t p1, const population_t p2) {
    std::list<reaction_t> rxn_list;
    if ((p1.species.name == "void") && (p2.species.name == "void")) {
        rxn_list.push_back(spawn_rxn(species_utilities::make_h_species()));
        rxn_list.push_back(spawn_rxn(species_utilities::make_p_species()));
    }
    return rxn_list;
}

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
        new_relation.reactions = get_reactions_for_species_pair(p1,p2);

        //make a relation_address_t held by the _second_ member of the pair
        relation_address_t new_relation_address;
        new_relation_address.owner_population_ptr = &p2;
        new_relation_address.relation_ptr = &new_relation;

        //put them in the right places
        p1.relations.push_back(new_relation);
        p2.relation_addresses.push_back(new_relation_address);
    }
}

ensemble_t initialize_ensemble() {
    //start at time zero with zero propensity
    ensemble_t ensemble;
    ensemble.current_time = 0.0;
    ensemble.total_propensity = 0.0;

    //ensembles always start with the void population
    add_population_to_ensemble(ensemble,pop_utilities::make_void_population());

    return ensemble;
}

int main() {
    auto ensemble = initialize_ensemble();
}

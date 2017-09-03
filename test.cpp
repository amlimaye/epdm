#include "types.hpp"
#include "constants.hpp"
#include "utilities.hpp"
#include <iostream>
#include <numeric>
#include <random>
#include <cmath>

void add_population_to_ensemble(ensemble_t& ensemble, population_t&& population) {
    //add the population to the ensemble
    ensemble.populations.push_back(population);

    //make tuples representing all _new_ bimolecular reactions with this species
    std::list<std::tuple<population_t&, population_t&>> all_pairs;
    for (auto& pop : ensemble.populations) {
        /*
            this caused me a few minutes of grief, so leaving a detailed comment.
            std::make_tuple makes tuples of values by default 
                (https://stackoverflow.com/questions/19054347)
            here, we actually want tuples of _references_ since we are going to
            fiddle with the members of the population_t struct inside the for loop.
            one forces creation of a tuple of _references_ with std::ref()
        */
        all_pairs.push_back(std::make_tuple(std::ref(pop),std::ref(population)));
    }

    //create relations for each pair
    for (auto& pair : all_pairs) {
        /*
            another few minutes of grief here.
            the std::get<>() in this function is out of control. AFAICT there isn't
            a way to unpack _references_ from a tuple in C++11. Structured bindings in
            C++17 solve this issue, but I chose not to use that to maintain C++11
            compatibility
        */

        //make a relation_t object held by the _first_ member of the pair
        relation_t new_relation;
        new_relation.owner_population_ptr = &std::get<0>(pair);

        //add reactions for this species pair
        new_relation.reactions = rxn_utilities::get_reactions_for_species_pair(std::get<0>(pair),std::get<1>(pair));
        new_relation.tot_partial_propensity = 0.0;

        //set the total partial propensity for this relation
        for (auto rxn : new_relation.reactions) {
            new_relation.tot_partial_propensity += rxn.partial_propensity;
        }

        //make a relation_address_t held by the _second_ member of the pair
        relation_address_t new_relation_address;
        new_relation_address.owner_population_ptr = &std::get<1>(pair);
        new_relation_address.relation_ptr = &new_relation;

        //put them in the right places
        std::get<0>(pair).relations.push_back(new_relation);
        std::get<1>(pair).relation_addresses.push_back(new_relation_address);

        //update the propensities of p1
        std::get<0>(pair).tot_partial_propensities.push_back(new_relation.tot_partial_propensity);
        std::get<0>(pair).tot_propensity = std::accumulate(std::get<0>(pair).tot_partial_propensities.begin(),
                                                           std::get<0>(pair).tot_partial_propensities.end(),
                                                           0.0);
        std::get<0>(pair).tot_full_propensity = std::get<0>(pair).tot_propensity*std::get<0>(pair).num_molecules;
    }

    //update the total propensity of the ensemble
    ensemble.total_propensity = 0.0;
    for (auto pop : ensemble.populations) {
        ensemble.total_propensity += pop.tot_full_propensity;
    }
}

void delete_population_from_ensemble(ensemble_t& ensemble) {
}

const reaction_t& sample_reaction(const ensemble_t& ensemble, double rand) {
    //select a population from the ensemble's list
    auto pop_it = ensemble.populations.begin();
    double s1 = 0;
    double s2 = 0;
    while (ensemble.total_propensity*rand > s1) {
        s2 = s1;
        s1 += pop_it->tot_full_propensity;
        std::advance(pop_it,1);
    }
    pop_it = std::prev(pop_it);

    //select a relation from the selected population's list    
    double g = (ensemble.total_propensity*rand - s2)/(pop_it->num_molecules);
    auto rel_it = pop_it->relations.begin();
    while (g > 0) {
        g -= rel_it->tot_partial_propensity;
        std::advance(rel_it,1);
    }
    rel_it = std::prev(rel_it);
    g += rel_it-> tot_partial_propensity;

    //select a reaction from this relation's list
    auto rxn_it = rel_it->reactions.begin();
    while (g > 0) {
        g -= rxn_it->partial_propensity;
        std::advance(rxn_it,1);
    }
    rxn_it = std::prev(rxn_it);

    //return a reference to the chosen reaction
    return *rxn_it;
}

void take_timestep(ensemble_t& ensemble) {
    //get two random numbers from the real uniform distribution over [0,1)
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    double r1 = distribution(ensemble.generator);
    double r2 = distribution(ensemble.generator);

    //sample the next time and the next reaction
    double next_time = (1.0/ensemble.total_propensity)*std::log((1.0/r1));
    auto next_rxn = sample_reaction(ensemble,r2);
    rxn_utilities::print_reaction(next_rxn);
}

ensemble_t initialize_ensemble() {
    //start at time zero with zero propensity
    ensemble_t ensemble;
    ensemble.current_time = 0.0;
    ensemble.total_propensity = 0.0;

    //ensembles always start with the void population
    add_population_to_ensemble(ensemble,pop_utilities::make_void_population());

    for (auto elem : ensemble.populations) {
        pop_utilities::print_population(elem);
    }

    return ensemble;
}

int main() {
    auto ensemble = initialize_ensemble();
    take_timestep(ensemble);
}

#include <iostream>
#include "../include/species_utilities.hpp"
#include "../include/pop_utilities.hpp"
#include "../include/rxn_utilities.hpp"

population_t pop_utilities::make_void_population() {
    //make the void species first
    auto void_species = species_utilities::make_void_species();

    //make a void population that owns the species
    population_t void_pop;
    void_pop.species = void_species;
    void_pop.num_molecules = 1;
    void_pop.tot_propensity = 0;

    return void_pop;
}

void pop_utilities::update_molecule_count(population_t* pop, long int delta) {
    //increment its molecule count
    pop->num_molecules += delta;

    //need to recompute the partial propensity in each relation that involves the incremented species
    for (auto& ra : pop->relation_addresses) {
        //first subtract off the current tot_partial propensity and set it to zero, we are recomputing it
        ra.relation_ptr->owner_population_ptr->tot_propensity -= ra.relation_ptr->tot_partial_propensity;
        ra.relation_ptr->tot_partial_propensity = 0;

        //recompute partial_propensity for each reaction wrt its owner
        for (auto& rxn : ra.relation_ptr->reactions) {
            double new_pp = rxn_utilities::compute_partial_propensity(rxn,pop->num_molecules);
            rxn.partial_propensity = new_pp;
            ra.relation_ptr->tot_partial_propensity += new_pp;
        }

        //add back the freshly recomputed tot_partial_propensity and also set the tot_full_propensity to the right value
        ra.relation_ptr->owner_population_ptr->tot_propensity += ra.relation_ptr->tot_partial_propensity;
        ra.relation_ptr->owner_population_ptr->tot_full_propensity = 
                ra.relation_ptr->owner_population_ptr->tot_propensity * ra.relation_ptr->owner_population_ptr->num_molecules;
    }
}

void pop_utilities::print_population(const population_t& in){
	std::cout << "Population {" << in.species.name << "}"<< std::endl;
	std::cout << "\tn = " << in.num_molecules << std::endl;
}
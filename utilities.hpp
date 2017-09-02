#include "types.hpp"
#include <iostream>
#include <tuple>

namespace species_utilities {
	species_t make_void_species() {
	    species_t void_species;
	    void_species.name = "void";
	    return void_species;
	}

	species_t make_h_species() {
	    species_t h_species;
	    h_species.name = "H";
	    return h_species;
	}

	species_t make_p_species() {
	    species_t p_species;
	    p_species.name = "P";
	    return p_species;
	}	
}

namespace pop_utilities {
	population_t make_void_population() {
	    //make the void species first
	    auto void_species = species_utilities::make_void_species();

	    //make a void population that owns the species
	    population_t void_pop;
	    void_pop.species = void_species;
	    void_pop.num_molecules = 1;

	    return void_pop;
	}

	void print_population(const population_t& in){
		std::cout << "Population {" << in.species.name << "}"<< std::endl;
		std::cout << "\tn = " << in.num_molecules << std::endl;
	}
}

namespace rxn_utilities {
	reaction_t spawn_rxn(const species_t spawned_species) {
	    reaction_t new_reaction;

	    //for a spawning reaction, the partial propensity and the rate constant are the same
	    new_reaction.rate_constant = constants::alpha;
	    new_reaction.partial_propensity = constants::alpha;

	    //add in products and reactants with stoichiometries
	    std::list<std::tuple<species_t, long int>> reactants, products;
	    reactants.push_back(std::make_tuple(species_utilities::make_void_species(),0));
	    reactants.push_back(std::make_tuple(species_utilities::make_void_species(),0));
	    products.push_back(std::make_tuple(spawned_species,1));
	    new_reaction.reactants = reactants;
	    new_reaction.products = products;

	    return new_reaction;
	}

	std::list<reaction_t> get_reactions_for_species_pair(const population_t p1, const population_t p2) {
	    std::list<reaction_t> rxn_list;
	    if ((p1.species.name == "void") && (p2.species.name == "void")) {
	        rxn_list.push_back(rxn_utilities::spawn_rxn(species_utilities::make_h_species()));
	        rxn_list.push_back(rxn_utilities::spawn_rxn(species_utilities::make_p_species()));
	    }
	    return rxn_list;
	}

	void print_reaction(const reaction_t& in) {
		std::cout << "Reactants:" << std::endl;
		for (auto rxt : in.reactants) {
			std::cout << std::get<1>(rxt) << " " << std::get<0>(rxt).name << std::endl;
		}
		std::cout << "Products:" << std::endl;
		for (auto prd : in.products) {
			std::cout << std::get<1>(prd) << " " << std::get<0>(prd).name << std::endl;
		}
	}
}
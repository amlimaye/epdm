#include "types.hpp"
#include <iostream>
#include <tuple>
#include <json/json.h>

namespace species_utilities {
	species_t make_void_species() {
	    species_t void_species;
	    void_species.name = "void";
        void_species.folded = false;
	    return void_species;
	}

	species_t make_h_species() {
	    species_t h_species;
	    h_species.name = "H";
        h_species.folded = false;
	    return h_species;
	}

	species_t make_p_species() {
	    species_t p_species;
	    p_species.name = "P";
        p_species.folded = false;
	    return p_species;
	}	

	species_t make_arbitrary_species(std::string name) {
		species_t new_species;
		new_species.name = name;
        new_species.folded = false;
		return new_species;
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
	    void_pop.tot_propensity = 0;

	    return void_pop;
	}

	void print_population(const population_t& in){
		std::cout << "Population {" << in.species.name << "}"<< std::endl;
		std::cout << "\tn = " << in.num_molecules << std::endl;
	}
}

namespace ensemble_utilities {
	void print_ensemble(const ensemble_t& in) {
		std::cout << "time = " << in.current_time << std::endl;
	    for (auto pop : in.populations) {
            if (pop.species.folded)
	    	    std::cout << "\t " << pop.species.name << "_f" << ": " << pop.num_molecules << std::endl;
            else
	    	    std::cout << "\t " << pop.species.name << "_u" << ": " << pop.num_molecules << std::endl;
	    }
	}

	Json::Value serialize_to_json(const ensemble_t& in) {
		Json::Value out;
		out["current_time"] = in.current_time;

		for (auto pop : in.populations) {
			out["populations"][pop.species.name] = Json::Int64(pop.num_molecules);
		}

		return out;
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

	reaction_t decay_rxn(const species_t decay_species) {
		reaction_t new_reaction;

		new_reaction.rate_constant = constants::delta;
	    new_reaction.partial_propensity = constants::delta;

	    //add in products and reactants with stoichiometries
	    std::list<std::tuple<species_t, long int>> reactants, products;
	    reactants.push_back(std::make_tuple(species_utilities::make_void_species(),0));
	    reactants.push_back(std::make_tuple(decay_species,1));
	    products.push_back(std::make_tuple(species_utilities::make_void_species(),0));
	    new_reaction.reactants = reactants;
	    new_reaction.products = products;

	    return new_reaction;
	}

	reaction_t elongation_rxn(const species_t s1, const species_t s2) {
		reaction_t new_reaction;
		new_reaction.rate_constant = constants::beta;

		//add in products and reactants with stoichiometries
	    std::list<std::tuple<species_t, long int>> reactants, products;
	    reactants.push_back(std::make_tuple(s1,1));
	    reactants.push_back(std::make_tuple(s2,1));
	    products.push_back(std::make_tuple(species_utilities::make_arbitrary_species(s1.name + s2.name),1));

	    new_reaction.reactants = reactants;
	    new_reaction.products = products;

	    return new_reaction;
	}

	reaction_t splitting_rxn(const species_t s1, size_t position) {
		reaction_t new_reaction;
		new_reaction.rate_constant = constants::hyd;

		//add in products and reactants with stoichiometries
	    std::list<std::tuple<species_t, long int>> reactants, products;
	    reactants.push_back(std::make_tuple(species_utilities::make_void_species(),0));
	    reactants.push_back(std::make_tuple(s1,1));
	    products.push_back(std::make_tuple(species_utilities::make_arbitrary_species(s1.name.substr(0,position+1)),1));
	    products.push_back(std::make_tuple(species_utilities::make_arbitrary_species(s1.name.substr(position+1)),1));

	    new_reaction.reactants = reactants;
	    new_reaction.products = products;

	    return new_reaction;
	}

    reaction_t unfolding_rxn(species_t in) {
        reaction_t new_reaction;
		new_reaction.rate_constant = constants::unfold;

		//add in products and reactants with stoichiometries
	    std::list<std::tuple<species_t, long int>> reactants, products;
	    reactants.push_back(std::make_tuple(species_utilities::make_void_species(),0));
	    reactants.push_back(std::make_tuple(in,1));

        in.folded = false;
        products.push_back(std::make_tuple(in,1));

	    new_reaction.reactants = reactants;
	    new_reaction.products = products;

	    return new_reaction;
    } 

    reaction_t folding_rxn(species_t in) {
        reaction_t new_reaction;
		new_reaction.rate_constant = constants::fold;

		//add in products and reactants with stoichiometries
	    std::list<std::tuple<species_t, long int>> reactants, products;
	    reactants.push_back(std::make_tuple(species_utilities::make_void_species(),0));
	    reactants.push_back(std::make_tuple(in,1));

        in.folded = true;
        products.push_back(std::make_tuple(in,1));

	    new_reaction.reactants = reactants;
	    new_reaction.products = products;

	    return new_reaction;
    } 

	double compute_partial_propensity(const reaction_t& reaction, long int num_molecules) {
		auto first_rxtnt_it = reaction.reactants.begin();
		auto second_rxtnt_it = std::next(first_rxtnt_it,1);

		auto first_rxtnt_name = std::get<0>(*first_rxtnt_it).name;
		auto second_rxtnt_name = std::get<0>(*second_rxtnt_it).name;

		if ((first_rxtnt_name == std::string("void")) || (second_rxtnt_name == std::string("void"))) {
			//this is either unimolecular or a source reaction. either way, pp = rate constant
			return reaction.rate_constant*num_molecules;
		} else if (first_rxtnt_name == second_rxtnt_name) {
			//homogenous bimolecular reaction
			return reaction.rate_constant*(num_molecules-1)*0.5;
		} else {
			//heterogeneous bimolecular reaction
			return reaction.rate_constant*num_molecules;
		}
	}

	std::list<reaction_t> get_reactions_for_species_pair(const population_t& p1, const population_t& p2) {
	    std::list<reaction_t> rxn_list;

	    // add: [V + V -> H], [V + V -> P]
	    if ((p1.species.name == "void") && (p2.species.name == "void")) {
	    	auto h_rxn = rxn_utilities::spawn_rxn(species_utilities::make_arbitrary_species("H"));
	    	auto p_rxn = rxn_utilities::spawn_rxn(species_utilities::make_arbitrary_species("P"));

	        rxn_list.push_back(h_rxn);
	        rxn_list.push_back(p_rxn);
	    } 

	    // add: [V + P -> V], [V + H -> V] 
	    // note, it turns out that the void species will always be the first one here
	    if ((p1.species.name == "void") && (p2.species.name == "P" || p2.species.name == "H")) {
	    	rxn_list.push_back(rxn_utilities::decay_rxn(species_utilities::make_arbitrary_species(p2.species.name)));
	    }

        // add: [V + {seq}_u -> {seq}_f]
        if ((p1.species.name == "void") && (p2.species.name != "void") && (p2.species.name.length() > 2) && (!p2.species.folded)) {
            rxn_list.push_back(rxn_utilities::folding_rxn(p2.species));
            #ifdef __DEBUG
            std::cout << "added a folding reaction for " << p2.species.name << std::endl;
            #endif
        }

        // add: [V + {seq}_u -> {seq}_f]
        if ((p1.species.name == "void") && (p2.species.name != "void") && (p2.species.folded)) {
            rxn_list.push_back(rxn_utilities::unfolding_rxn(p2.species));
            #ifdef __DEBUG
            std::cout << "added an unfolding reaction for " << p2.species.name << std::endl;
            #endif
        }

	    // add: polymerization reaction
	    if ((p1.species.name.length() == 1) && (p2.species.name.length() >= 1) && (p2.species.name != "void") && (!p2.species.folded)) {
	    	rxn_list.push_back(rxn_utilities::elongation_rxn(species_utilities::make_arbitrary_species(p1.species.name),species_utilities::make_arbitrary_species(p2.species.name)));
	    	rxn_list.push_back(rxn_utilities::elongation_rxn(species_utilities::make_arbitrary_species(p2.species.name),species_utilities::make_arbitrary_species(p1.species.name)));
	    }

	    // add: polymerization reaction
		if ((p2.species.name.length() == 1) && (!p1.species.folded) && (p1.species.name.length() >= 1) && (p1.species.name != "void")) {
	    	rxn_list.push_back(rxn_utilities::elongation_rxn(species_utilities::make_arbitrary_species(p1.species.name),species_utilities::make_arbitrary_species(p2.species.name)));
	    	rxn_list.push_back(rxn_utilities::elongation_rxn(species_utilities::make_arbitrary_species(p2.species.name),species_utilities::make_arbitrary_species(p1.species.name)));
	    }

	    // add: splitting reaction(s)
	    if ((p1.species.name == "void") && ((p2.species.name != "void") && (p2.species.name.length() > 1))) {
	    	for (int pos = 0; pos < p2.species.name.length()-1; pos++) {
	    		rxn_list.push_back(rxn_utilities::splitting_rxn(p2.species,pos));
	    	}
	    	#ifdef __DEBUG
	    	std::cout << "added " << p2.species.name.length()-1 << " splitting reactions!"<< std::endl;
	    	#endif
	    }

	    return rxn_list;
	}

	void print_reaction(const reaction_t& in) {
		auto it = in.reactants.begin();
		for (int i = 0; i < in.reactants.size()-1; i++) {
            if (std::get<0>(*it).folded)
			    std::cout << std::get<1>(*it) << " " << std::get<0>(*it).name << "_f" << " + ";
            else
			    std::cout << std::get<1>(*it) << " " << std::get<0>(*it).name << "_u" << " + ";
			std::advance(it,1);
		}
        
        if (std::get<0>(*it).folded)
		    std::cout << std::get<1>(*it) << " " << std::get<0>(*it).name << "_f" << " --> ";
        else
		    std::cout << std::get<1>(*it) << " " << std::get<0>(*it).name << "_u" << " --> ";

		auto it2 = in.products.begin();
		for (int i = 0; i < in.products.size()-1; i++) {
            if (std::get<0>(*it2).folded)
			    std::cout << std::get<1>(*it2) << " " << std::get<0>(*it2).name << "_f" << " + ";
            else
			    std::cout << std::get<1>(*it2) << " " << std::get<0>(*it2).name << "_u" << " + ";
			std::advance(it2,1);
		}

        if (std::get<0>(*it2).folded)
		    std::cout << std::get<1>(*it2) << " " << std::get<0>(*it2).name << "_f" << std::endl;
        else
		    std::cout << std::get<1>(*it2) << " " << std::get<0>(*it2).name << "_u" << std::endl;
	}
}

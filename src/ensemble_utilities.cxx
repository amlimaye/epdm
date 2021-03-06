#include <iostream>
#include <random>
#include "pop_utilities.hxx"
#include "rxn_utilities.hxx"
#include "ensemble_utilities.hxx"
#include "folding.hxx"

ensemble_t ensemble_utilities::initialize_ensemble(const uint32_t seed) {
    //start at time zero with zero propensity
    ensemble_t ensemble;
    ensemble.current_time = 0.0;
    ensemble.total_propensity = 0.0;
    ensemble.generator.seed(seed);

    //ensembles always start with the void population
    add_population_to_ensemble(&ensemble,pop_utilities::make_void_population());

    return ensemble;	
}

void ensemble_utilities::add_population_to_ensemble(ensemble_t* ensemble, population_t population) {
    #ifdef __INCLUDE_FOLDING
    //determine whether or not this sequence is foldable
    bool foldable = false;
    int contacts = 0;
    int max_count = 0;
    if (!(population.species.name.length() > 13))
        std::tie(foldable,contacts,max_count) = folding::is_foldable(population.species.name);

    if (foldable) {
        #ifdef __DEBUG
        std::cout << "sequence " << population.species.name << " is foldable with " << contacts << " hydrophobic contacts, multiplicity " << max_count << "!" << std::endl;
        #endif
        population.species.foldable = true;
        population.species.native_contacts = contacts;
    } else {
        #ifdef __DEBUG
        std::cout << "sequence " << population.species.name << " is not foldable with " << contacts << " hydrophobic contacts, multiplicity " << max_count << "!" << std::endl;
        #endif
    }
    #endif

    //add the population to the ensemble
    ensemble->populations.push_back(population);

    //make tuples representing all _new_ bimolecular reactions with this species
    std::list<std::tuple<population_t*, population_t*>> all_pairs;
    for (auto& pop : ensemble->populations) {
        all_pairs.push_back(std::forward_as_tuple(&pop,&(ensemble->populations.back())));
    }

    //create relations for each pair
    for (auto& pair : all_pairs) {
        population_t* first_elem = std::get<0>(pair);
        population_t* second_elem = std::get<1>(pair);

        //make a relation_t object held by the _first_ member of the pair
        relation_t new_relation;
        new_relation.owner_population_ptr = first_elem;

        //add reactions for this species pair
        new_relation.reactions = rxn_utilities::get_reactions_for_species_pair(*first_elem,*second_elem);
        if (new_relation.reactions.size() != 0) {
            new_relation.tot_partial_propensity = 0.0;

            //set the total partial propensity for this relation
            for (auto rxn : new_relation.reactions) {
                new_relation.tot_partial_propensity += rxn.partial_propensity;
            }

            //put them in the right places
            first_elem->relations.push_back(new_relation);

            //make a relation_address_t held by the _second_ member of the pair
            relation_address_t new_relation_address;
            new_relation_address.owner_population_ptr = second_elem;
            new_relation_address.relation_ptr = &(first_elem->relations.back());

            second_elem->relation_addresses.push_back(new_relation_address);
            first_elem->relations.back().relation_address_ptr = &(second_elem->relation_addresses.back());

            //update the propensities of p1
            first_elem->tot_propensity += new_relation.tot_partial_propensity;
            first_elem->tot_full_propensity = (first_elem->tot_propensity)*(first_elem->num_molecules);
        }/* else {
            std::cout << "~~ added no reactions for pair: " << first_elem->species.name << ", " << second_elem->species.name << std::endl;
        }*/
    }

    //update the total propensity of the ensemble
    ensemble->total_propensity = 0.0;
    for (auto pop : ensemble->populations) {
        ensemble->total_propensity += pop.tot_full_propensity;
    }
}

const reaction_t ensemble_utilities::sample_reaction(const ensemble_t& ensemble, double rand) {
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

    //return the chosen reaction
    return *rxn_it;
}

std::tuple<bool,size_t> 
ensemble_utilities::is_in_tabulated_populations(const ensemble_t& ensemble, 
												const species_t& species) {
    size_t position = 0;
    for (auto pop : ensemble.populations) {
        if ((species.name == pop.species.name) && (species.folded == pop.species.folded))
            return std::make_tuple(true,position);
        position++;
    }
    return std::make_tuple(false,position);
}

void ensemble_utilities::take_timestep(ensemble_t& ensemble) {
    //get two random numbers from the real uniform distribution over [0,1)
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    double r1 = distribution(ensemble.generator);
    double r2 = distribution(ensemble.generator);

    //sample the next time and the next reaction
    double next_time = (1.0/ensemble.total_propensity)*std::log((1.0/r1));
    ensemble.current_time += next_time;
    auto next_rxn = sample_reaction(ensemble,r2);

    #ifdef __DEBUG
    rxn_utilities::print_reaction(next_rxn);
    std::cout << ensemble.total_propensity << std::endl;
    #endif

    //loop through each of the products
    for (auto prod_tup : next_rxn.products) {
        auto prod_species = std::get<0>(prod_tup);
        auto prod_stoich = std::get<1>(prod_tup);

        size_t position;
        bool found;
        std::tie(found,position) = ensemble_utilities::is_in_tabulated_populations(ensemble,prod_species);

        //if this species is already being tabulated, we need to change its num_molecules
        //and then recompute partial propensities in each of the relations pointed to by its
        //relation addresses
        if (found) {
            //find the population in the ensemble
            auto pop_it = ensemble.populations.begin();
            std::advance(pop_it,position);

            pop_utilities::update_molecule_count(&(*pop_it),prod_stoich);
            //                    ^^^ ugly way to get a pointer from an iterator
        }
        //if this is a species we are _not_ yet tabulating, add the population to the ensemble
        else {
            population_t new_pop;
            new_pop.species = prod_species;
            new_pop.num_molecules = prod_stoich;
            new_pop.tot_propensity = 0;
            add_population_to_ensemble(&ensemble,new_pop);
        }
    }
    for (auto rxtnt_tup : next_rxn.reactants) {
        auto rxtnt_species = std::get<0>(rxtnt_tup);
        auto rxtnt_stoich = std::get<1>(rxtnt_tup);

        //the reactant is guaranteed to be in the tabulated population, but we still need its position
        size_t position;
        bool found;
        std::tie(found,position) = ensemble_utilities::is_in_tabulated_populations(ensemble,rxtnt_species);
        auto pop_it = ensemble.populations.begin();
        std::advance(pop_it,position);
        pop_utilities::update_molecule_count(&(*pop_it),-1*rxtnt_stoich);
        //                    ^^^ ugly way to get a pointer from an iterator

        //if this causes the molecule count to go to zero, we need to eliminate the population from the ensemble
        if (pop_it->num_molecules == 0) {
            //for all relations, delete any relation addresses that point to them
            for (auto& rel : pop_it->relations) {
                rel.relation_address_ptr->owner_population_ptr->relation_addresses.remove(*(rel.relation_address_ptr));
            }

            //for all relation addresses, delete the relation they point to
            for (auto& rel_addr : pop_it->relation_addresses) {
                rel_addr.relation_ptr->owner_population_ptr->relations.remove(*(rel_addr.relation_ptr));
            }

            //remove this population from the ensemble's list
            ensemble.populations.erase(pop_it);
        }
    }

    //update the total propensity of the ensemble
    ensemble.total_propensity = 0.0;
    for (auto& pop : ensemble.populations) {
        pop.tot_propensity = 0.0;
        for (auto relation : pop.relations) {
            pop.tot_propensity += relation.tot_partial_propensity;
        }
        ensemble.total_propensity += pop.tot_full_propensity;
    }

    #ifdef __DEBUG
    std::cout << ensemble.total_propensity << std::endl;
    #endif
}

void ensemble_utilities::print_ensemble(const ensemble_t& in) {
	std::cout << "time = " << in.current_time << std::endl;
    for (auto pop : in.populations) {
        if (pop.species.folded)
    	    std::cout << "\t " << pop.species.name << "_f" << ": " << pop.num_molecules 
    		<< std::endl;
        else
    	    std::cout << "\t " << pop.species.name << "_u" << ": " << pop.num_molecules 
    		<< std::endl;
    }
}

std::ostream& operator<<(std::ostream& stream, const ensemble_t& in) {
    stream << "time = " << in.current_time << std::endl;
    for (auto pop : in.populations) {
        if (pop.species.folded)
            stream << "\t " << pop.species.name << "_f" << ": " << pop.num_molecules 
            << std::endl;
        else
            stream << "\t " << pop.species.name << "_u" << ": " << pop.num_molecules 
            << std::endl;
    } 
    return stream;
}

Json::Value ensemble_utilities::serialize_to_json(const ensemble_t& in) {
	Json::Value out;
	out["current_time"] = in.current_time;

	for (auto pop : in.populations) {
		if (pop.species.folded)
			out["populations"][pop.species.name + "_f"] = Json::Int64(pop.num_molecules);
		else
			out["populations"][pop.species.name + "_u"] = Json::Int64(pop.num_molecules);
	}

	return out;
}

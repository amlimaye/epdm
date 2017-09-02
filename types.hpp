#ifndef __TYPES_H
#define __TYPES_H

#include <string>
#include <list>
#include <tuple>
#include <random>

//need to forward-declare some structs to avoid circular dependencies
struct population_t;
struct relation_address_t;

typedef struct species_t {
    std::string name;
} species_t;

typedef struct reaction_t {
    double                                      rate_constant;          //rate constant k for bimolecular rxn.
    double                                      partial_propensity;     //rate = k*n_1*n_2; propensity = rate/n_{owning_species}
    std::list<std::tuple<species_t, long int>>  reactants;              //[(species, stoich. coeff.), ...]
    std::list<std::tuple<species_t, long int>>  products;               //^same as above
} reaction_t;

typedef struct relation_t {
    std::list<reaction_t>   reactions;                                  //list of all possible reactions between these two species
    double                  tot_partial_propensity;                     //sum reaction.partial_propensity over all reactions
    population_t*           owner_population_ptr;                       //pointer to the population that "owns" this relation
    relation_address_t*     relation_address_ptr;                       //pointer to the relation address pointing to this relation, which will be held by the other species in this relation
} relation_t;

typedef struct relation_address_t {
    population_t*           owner_population_ptr;                       //pointer to the population that owns this relation address
    relation_t*             relation_ptr;                               //pointer to a relation
} relation_address_t;

typedef struct population_t {
    species_t                       species;                            //"owner species" of this population
    long int                        num_molecules;                      //total number of molecules of this population currently present
    double                          tot_propensity;                     //sum_i pp_i for all relations owned by this population
    std::list<double>               tot_partial_propensities;           //[pp_1, ...] for all relations owned by this population
    std::list<relation_t>           relations;                          //list of relations owned by this population
    std::list<relation_address_t>   relation_addresses;                 //list of relation_addresses owned by other populations, but pointing to this population
} population_t;

typedef struct ensemble_t {
    double                      current_time;
    double                      total_propensity;
    std::list<population_t>     populations;
    std::default_random_engine  generator;
} ensemble_t;

#endif
#include "../include/types.hpp"
#include "../include/constants.hpp"
#include "../include/species_utilities.hpp"
#include "../include/pop_utilities.hpp"
#include "../include/rxn_utilities.hpp"
#include "../include/ensemble_utilities.hpp"

#include <iostream>
#include <numeric>
#include <random>
#include <cmath>
#include <fstream>
#include <json/json.h>

Json::Value initialize_json() {
    Json::Value root;
    Json::Value time_list(Json::arrayValue);
    root["epdm"] = time_list;
    return root;
}

void dump_json(const Json::Value& json, const std::string& fname) {
    std::ofstream fout(fname);
    fout << json;
}

int main(int argc, char *argv[]) {
    //argument parsing...
    long int num_steps;
    uint32_t rand_seed;
    if (argc == 3) {
        num_steps = atol(argv[1]);
        rand_seed = atoi(argv[2]);
    } else {
        std::cout << "usage: " << argv[0] << " <num_steps> <rand_seed>" << std::endl;
        exit(1);
    }

    auto json = initialize_json();
    auto ensemble = ensemble_utilities::initialize_ensemble(rand_seed);

    //ensemble_utilities::print_ensemble(ensemble);
    for (int i = 0; i < num_steps; i++) {
        ensemble_utilities::take_timestep(ensemble);
        json["epdm"].append(ensemble_utilities::serialize_to_json(ensemble));

        #ifdef __DEBUG
        ensemble_utilities::print_ensemble(ensemble);
        #endif

        if ((i+1) % 100 == 0) {
            std::cout << "step " << i+1 << "/" << num_steps << std::endl;
        }
    }

    dump_json(json,"out.json");
}

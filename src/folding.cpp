#include "folding.hpp"
#include <iostream>
#include <fstream>
#include <json/json.h>

std::tuple<bool,int,int> folding::is_foldable(const std::string name) {
    //check some dumb things
    if (name.length() < 4)
        return std::make_tuple(false,0,0);
    if (name == "void")
        return std::make_tuple(false,0,0);

    //load up the contacts JSON file
    std::string fname = "contacts/contacts_" + std::to_string(name.length()) + ".json";
    #ifdef __DEBUG
    std::cout << fname << std::endl;
    #endif
    std::ifstream fin(fname);
    Json::Value root;
    fin >> root;

    //find the hydrophobic indices in _this_ sequence
    std::vector<int> hyd_indices;
    for (size_t i = 0; i < name.length(); i++) {
        if (name.at(i) == 'H')
            hyd_indices.push_back(i);
    }

    //bookkeeping variables
    int max_count = 0;
    int max_contacts = 0;

    //iterate through each configuration in the JSON
    for (auto& elem : root) {
        int contacts = 0;

        //iterate through the hydrophobic indices in each sequence
        for (auto hyd_idx : hyd_indices) {
            auto this_elem = elem[hyd_idx];

            //iterate through the neighbors of this index at which there is an H residue in the reference sequence
            for (auto neigh_idx : this_elem) {

                //add a contact if this sequence has an H residue at that neighbor position
                if (name.at(neigh_idx.asInt()) == 'H')
                    contacts += 1;
            }
        }

        //if we found a new max, set it and set its count to zero
        if (contacts > max_contacts) {
            max_count = 0;
            max_contacts = contacts;
        }

        //if we found a sequence with just as many contacts up the counter
        if (contacts == max_contacts) {
            max_count += 1;
        }
    }

    //foldable if it has at least one contact and has a unique native structure
    if ((max_count == 1) && (max_contacts > 1))
        return std::make_tuple(true,max_contacts,max_count);

    //else it's not foldable
    return std::make_tuple(false,max_contacts,max_count);
}

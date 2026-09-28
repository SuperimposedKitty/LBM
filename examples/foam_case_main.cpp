#include "lbm/foam_case.hpp"

#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char **argv) {
    try {
        if (argc != 3 || std::string(argv[1]) != "-case")
            throw std::invalid_argument("Usage: " LBM_CASE_TOOL " -case <case-directory>");
        const std::filesystem::path path = argv[2];
        const std::string tool = LBM_CASE_TOOL;
        if (tool == "lbm_blockMesh")
            lbm::block_mesh_case(path);
        else if (tool == "lbm_setFields")
            lbm::set_fields_case(path);
        else
            lbm::solve_foam_case(path);
        return 0;
    } catch (const std::exception &e) {
        std::cerr << LBM_CASE_TOOL << ": " << e.what() << '\n';
        return 1;
    }
}

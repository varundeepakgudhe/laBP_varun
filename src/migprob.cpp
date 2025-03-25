#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include "json.hpp" // JSON library
#include "migprob.h"

using json = nlohmann::json;

// Function to read migration matrix from JSON
Matrix readMigrationMatrix(const std::string &filename, const std::string &selected_time)
{
    // Open the JSON file
    std::ifstream file(filename);
    if (!file)
    {
        std::cerr << "Error: Unable to open JSON file.\n";
        return {};
    }

    // Parse JSON
    json j;
    file >> j;

    // Map to store migration matrices with time keys
    std::map<std::string, Matrix> migration_data;

    // Extract matrices and store them in the map
    for (auto &[time, matrix] : j.items())
    {
        migration_data[time] = matrix.get<Matrix>();
    }

    // Check if the selected time exists in the data
    if (migration_data.find(selected_time) == migration_data.end())
    {
        std::cerr << "Error: Time key " << selected_time << " not found in JSON.\n";
        return {};
    }

    // Return the migration matrix for the selected time
    return migration_data[selected_time];
}

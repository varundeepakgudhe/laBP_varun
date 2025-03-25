#ifndef MIGPROB_H
#define MIGPROB_H

#include <vector>
#include <string>

// Define the Matrix type
using Matrix = std::vector<std::vector<double>>;

// Function to read the migration matrix from JSON
Matrix readMigrationMatrix(const std::string &filename, const std::string &selected_time);

#endif // MIGPROB_H

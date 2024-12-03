#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <iomanip>

// Function to trim whitespace from both ends of a string
std::string trim(const std::string& s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    size_t end = s.find_last_not_of(" \t\r\n");
    return (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
}

// Function to split a string by a given delimiter
std::vector<std::string> split(const std::string &s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(trim(token));
    }
    return tokens;
}

void TabulateActiveTriggerInformation() {
    // Path to your CSV file
    std::string filename = "/Users/patsfan753/Desktop/triggerAnalysisCombined_GoldenRunList_12_4_24.csv";
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "\033[31mError opening file: " << filename << "\033[0m" << std::endl;
        return;
    }

    std::string line;
    // Read the header line to get column names
    std::getline(file, line);

    // Determine the delimiter (comma or tab)
    std::vector<std::string> headers = split(line, ',');
    if (headers.size() <= 1) { // If only one column, try tab delimiter
        headers = split(line, '\t');
    }

    // Check if headers were properly read
    if (headers.size() <= 1) {
        std::cerr << "\033[31mError: Unable to parse headers. Please check the delimiter.\033[0m" << std::endl;
        return;
    }

    // Initialize counts for each trigger
    std::map<std::string, int> offCounts;
    std::map<std::string, int> onCounts;
    std::vector<std::string> triggerNames;

    for (size_t i = 1; i < headers.size(); ++i) { // Skip 'runNumber'
        offCounts[headers[i]] = 0;
        onCounts[headers[i]] = 0;
        triggerNames.push_back(headers[i]);
    }

    int totalRuns = 0;

    // Process each row in the CSV file
    while (std::getline(file, line)) {
        // Skip empty lines
        if (line.empty()) continue;

        totalRuns++;
        std::vector<std::string> tokens = split(line, ',');

        // If the number of tokens doesn't match, try splitting by tab
        if (tokens.size() != headers.size()) {
            tokens = split(line, '\t');
        }

        // Ensure data integrity
        if (tokens.size() != headers.size()) {
            std::cerr << "\033[33mWarning: Data format error in line " << totalRuns + 1 << " (expected "
                      << headers.size() << " tokens, got " << tokens.size() << ")\033[0m" << std::endl;
            continue;
        }

        // Count 'ON' and 'OFF' statuses for each trigger
        for (size_t i = 1; i < headers.size(); ++i) { // Skip 'runNumber'
            std::string value = tokens[i];
            if (value == "OFF") {
                offCounts[headers[i]]++;
            } else if (value == "ON") {
                onCounts[headers[i]]++;
            }
        }
    }

    // Output the results in a formatted table with ANSI color codes
    std::cout << "\n\033[1mTotal number of runs: " << totalRuns << "\033[0m\n" << std::endl;

    std::cout << "\033[1mTrigger Status Counts:\033[0m" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "| \033[1m" << std::setw(35) << std::left << "Trigger Name"
              << "\033[0m | \033[1m" << std::setw(10) << "OFF Count"
              << "\033[0m | \033[1m" << std::setw(10) << "ON Count" << "\033[0m |" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;

    for (const auto& trigger : triggerNames) {
        std::string offColor = (offCounts[trigger] == 0) ? "\033[32m" : "\033[31m"; // Green if zero, red otherwise
        std::string onColor = (onCounts[trigger] > 0) ? "\033[32m" : "\033[31m";    // Green if >0, red if zero

        std::cout << "| " << "\033[1m" << std::setw(35) << std::left << trigger << "\033[0m"
                  << " | " << offColor << std::setw(10) << offCounts[trigger] << "\033[0m"
                  << " | " << onColor << std::setw(10) << onCounts[trigger] << "\033[0m |" << std::endl;
    }

    std::cout << "--------------------------------------------------------------------" << std::endl;
}

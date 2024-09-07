#include <iostream>
#include <fstream>
#include <set>
#include <string>

void FinalizeGoldenRunList_calo_gl1_fem() {
    std::set<int> runs_file1, runs_file2, runs_file3, golden_runs;

    // Open the first file and read the run numbers into a set
    std::ifstream file1("FileLists/GoldenCalorimeterRunList.txt");
    if (!file1.is_open()) {
        std::cerr << "Error opening file: list_runnumber_calo.txt" << std::endl;
    }

    int run_number;
    while (file1 >> run_number) {
        runs_file1.insert(run_number);
    }
    file1.close();

    // Open the second file and read the run numbers into a set
    std::ifstream file2("FileLists/GoldenGL1RunList.txt");
    if (!file2.is_open()) {
        std::cerr << "Error opening file: GoldenGL1RunList.txt" << std::endl;
    }

    while (file2 >> run_number) {
        runs_file2.insert(run_number);
    }
    file2.close();

    // Open the third file and read the run numbers into a set
    std::ifstream file3("FileLists/GoldenFEMrunList.txt");
    if (!file3.is_open()) {
        std::cerr << "Error opening file: GoldenFEMrunList.txt" << std::endl;
    }

    while (file3 >> run_number) {
        runs_file3.insert(run_number);
    }
    file3.close();

    // Find golden run numbers (intersection of both sets)
    for (auto num : runs_file1) {
        if (runs_file2.find(num) != runs_file2.end() && runs_file3.find(num) != runs_file3.end()) {
            golden_runs.insert(num);
        }
    }

    // Save golden run numbers to a new file
    std::ofstream output_file("FileLists/FinalGoldenRunList_calo_gl1_femChecked.txt");
    if (!output_file.is_open()) {
        std::cerr << "Error opening file: FinalGoldenRunList_calo_gl1_femChecked.txt" << std::endl;
    }

    for (auto num : golden_runs) {
        output_file << num << std::endl;
    }
    output_file.close();

    std::cout << "Number of Calo + GL1 GOLDEN runs saved to FinalGoldenRunList_calo_gl1_femChecked.txt: " << golden_runs.size() << std::endl;
}

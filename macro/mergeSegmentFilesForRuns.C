#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>
#include <fstream>
#include <dirent.h>
#include <sys/stat.h>
#include <algorithm>
#include <regex>

// ROOT headers
#include <TSystem.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TObject.h>
#include <TKey.h>
#include <fcntl.h>    // For open
#include <unistd.h>   // For write, close

// Define ANSI color codes
#define ANSI_RESET  "\x1b[0m"
#define ANSI_RED    "\x1b[31m"
#define ANSI_GREEN  "\x1b[32m"
#define ANSI_YELLOW "\x1b[33m"
#define ANSI_CYAN   "\x1b[36m"

void get_scaledowns(int runnumber, int scaledowns[]) {
    // Connect to the PostgreSQL server
    TSQLServer *db = TSQLServer::Connect("pgsql://sphnxdaqdbreplica:5432/daq", "phnxro", "");

    // Check if the connection was successful
    if (db) {
        printf("Server info: %s\n", db->ServerInfo()); // Print server information
    } else {
        printf("Failed to connect to the database\n"); // Print error message if connection fails
        return; // Exit the function if the connection fails
    }

    // Initialize scaledowns to -1 (or any invalid value for your context)
    for (int i = 0; i < 64; i++) {
        scaledowns[i] = -1;
    }

    // Query to get all scaledown factors for the runnumber
    char sql[1000];
    sprintf(sql, "SELECT * FROM gl1_scaledown WHERE runnumber = %d;", runnumber);
    printf("Executing query: %s\n", sql); // Print the SQL query for debugging

    // Execute the query
    TSQLResult* res = db->Query(sql);

    // Check the number of rows returned
    int nrows = res->GetRowCount();
    if (nrows > 0) {
        TSQLRow* row = res->Next(); // There should be only one row for the runnumber
        int nfields = res->GetFieldCount();

        // Map column names to indices
        std::map<std::string, int> columnIndices;
        for (int i = 0; i < nfields; i++) {
            columnIndices[res->GetFieldName(i)] = i;
        }

        // For each bit, get the scaledown factor
        for (int is = 0; is < 64; is++) {
            char colName[20];
            sprintf(colName, "scaledown%02d", is);

            auto it = columnIndices.find(colName);
            if (it != columnIndices.end()) {
                int colIndex = it->second;
                const char* field = row->GetField(colIndex);
                if (field) {
                    scaledowns[is] = atoi(field);
                } else {
                    scaledowns[is] = -1;
                }
            } else {
                // Column not found
                scaledowns[is] = -1;
            }
        }

        delete row; // Clean up the row object
    } else {
        printf("No rows returned for runnumber %d\n", runnumber);
    }

    delete res; // Clean up the result set object
    delete db; // Close the database connection

    // Print out the scaledown factors for debugging
    printf("Scaledown factors for runnumber %d:\n", runnumber);
    for (int is = 0; is < 64; is++) {
        printf("  Scaledown factor for bit %d: %d\n", is, scaledowns[is]);
    }
}


void get_trigger_name_to_index_map(int runNumber, std::map<std::string, int>& triggerNameToIndexMap) {
    // Map1
    std::vector<int>  runNumbersForMap1 = { 44477,44478,44482,44483,44495,44498,44499,44503,44505,44506,44507,44509,44510,44511,44512,44513,44533,44534,44604,44608,44611,
        44616,44618,44619,44621,44631,44638,44642,45034,45035,45036,45038,45041,45048,45051,45052,45090,45100,45103,45105,45106,45107,
        45150,45151,45153,45154,45155,45157,45159,45160,45161,45162,45164,45166,45167,45170,45172,45176,45177,45178,45181,45183,45186,
        45189,45190,45191,45196,45199,45201,45203,45246,45248,45249,45252,45255,45256,45258,45274,45288,45290,45291,45292,45315,45316,
        45318,45325,45390,45391,45393,45394,45401,45402,45414,45443,45485,45486,45487,45489,45490,45491,45493,45494,45495,45531,45540,
        45541,45547,45548,45551,45552,45620,45624,45627,45628,45633,45637,45645,45724,45807,45816,45837,45841,45842,45851,45852,45856,
        45858,45872,45883,46011,46019,46022,46023,46025,46029,46036};
    

    std::map<int, std::string> triggerNameMap1 = {
        {0, "Clock"},
        {1, "ZDC_South"},
        {2, "ZDC_North"},
        {3, "ZDC_Coincidence"},
        {4, "HCAL_Singles"},
        {5, "HCAL_Coincidence"},
        {8, "MBD_S_geq_1"},
        {9, "MBD_N_geq_1"},
        {10, "MBD_NandS_geq_1"},
        {11, "MBD_NandS_geq_2"},
        {12, "MBD_NandS_geq_1_vtx_lessTehn_10_cm"},
        {13, "MBD_NandS_geq_1_vtx_lessThen_30_cm"},
        {14, "MBD_NandS_geq_1_vtx_lessThen_60_cm"},
        {15, "HCAL_Singles_plus_MBD_NS_geq_1"},
        {16, "Jet_4_GeV_plus_MBD_NS_geq_1"},
        {17, "Jet_6_GeV_plus_MBD_NS_geq_1"},
        {18, "Jet_8_GeV_plus_MBD_NS_geq_1"},
        {19, "Jet_10_GeV_plus_MBD_NS_geq_1"},
        {20, "Jet_4_GeV"},
        {21, "Jet_6_GeV"},
        {22, "Jet_8_GeV"},
        {23, "Jet_10_GeV"},
        {24, "Photon_1_GeV_plus_MBD_NS_geq_1"},
        {25, "Photon_2_GeV_plus_MBD_NS_geq_1"},
        {26, "Photon_3_GeV_plus_MBD_NS_geq_1"},
        {27, "Photon_4_GeV_plus_MBD_NS_geq_1"},
        {28, "Photon_1_GeV"},
        {29, "Photon_2_GeV"},
        {30, "Photon_3_GeV"},
        {31, "Photon_4_GeV"}
    };

    std::vector<int> runNumbersForMap2 = {
       46065,46068,46105,46107,46133,46135,46137,46427,46429,46430,46431,46432,46434,46436,46437,46438,46451,46452,46455,46473,46477,
       46480,46523,46524,46529,46530,46531,46537,46539,46541,46545,46546,46553,46559,46563,46567,46569,46577,46587,46588,46593,46595,
       46598,46603,46605,46619,46623,46640,46649,46655,46656,46663,46665,46668,46669,46676,46697,46699,46704,46721,46735,46750,46753,
       46755,46756,46759,46772,46775,46912,46917,46918,46941,46943,46944,46945,46946,46947,46949,46950,46951,46952,46965,46967,46968,
       47002,47006,47007,47009,47014,47017,47019,47022,47032,47033,47034,47036,47038,47040,47043,47051,47052,47053,47055,47056,47058,
       47060,47061,47064,47066,47068,47089,47098,47101,47102,47114,47115,47116,47124,47125,47129,47131,47135,47137,47138,47139,47140,
       47141,47143,47146,47155,47156,47158,47160,47161,47162,47201,47202,47203,47204,47211,47216,47219,47229,47230,47289,47293,47303,
       47306,47310,47315,47316,47323,47330,47332,47334,47360,47375,47376,47377,47378,47381,47382,47391,47393,47395,47396,47399,47443,
       47451,47455,47457,47458,47459,47464,47474,47476,47480,47484,47485,47491,47492,47494,47495,47497,47502,47503,47505,47506,47507,
       47513,47514,47516,47522,47524,47525,47538,47540,47548,47552,47557,47568,47634,47636,47638,47657,47658,47659,47661,47662,47666,
       47667,47698,47715,47716,47720,47722,47723,47724,47725,47727,47729,47730,47732,47733,47769,47777,47778,47783,47807,47831,47846,
       47848,47867,47893,47939,47946,47962,47966,47982,48073,48080,48085,48086,48100,48166,48180,48181,48231,48233,48234,48237,48239,
       48240,48244,48245,48253,48255,48256,48257,48258,48260,48261,48262,48263,48265,48287,48291,48293,48294,48295,48307,48313,48318,
       48320,48323,48325,48326,48327,48335,48337,48338,48341,48342,48343,48346,48347,48348,48349,48352,48356,48357,48358,48359,48366,
       48367,48369,48409,48410,48412,48416,48417,48418,48421,48422,48423,48454,48455,48456,48459,48461,48462,48469,48536,48635,48636,
       48638,48645,48656,48657,48658,48660,48701,48720,48721,48722,48725,48726,48727,48730,48731,48732,48734,48736,48742,48743,48745,
       48746,48801,48803,48805,48806,48807,48810,48811,48813,48824,48826,48828,48829,48832,48836,48838,48839,48859,48861,48863,48864,
       48865,48867,48868,48869,48870,48872,48873,48874,48877,48883,48884,48885,48891,48892,48893,48894,48895,48896,48897,48899,48900,
       48901,48902,48903,48918,48935,48936,48938,48940,48943,48946,48949,48951,48971,48976,48981,48982,48983,48984,48985,48986,48987,
       48988,48990,48991,49017,49023,49026,49027,49028,49029,49030,49031,49035,49042,49044,49047,49048,49050,49052,49053,49054,49060,
       49061,49062,49063,49066,49067,49069,49070,49071,49072,49073,49098,49125,49128,49133,49138,49219,49224,49226,49227,49228,49229,
       49230,49233,49240,49241,49244,49247,49248,49249,49250,49251,49254,49263,49264,49266,49267,49268,49269,49270,49307,49308,49309,
       49310,49311,49312,49313,49314,49316,49317,49324,49329,49330,49331,49332,49333,49336,49337,49338,49339,49340,49343,49345,49346,
       49347,49348,49349,49350,49351,49352,49356,49357,49358,49359,49361,49362,49363,49365,49366,49367,49368,49372,49374,49376,49377,
       49378,49379,49380,49381,49382,49383,49384,49385,49386,49389,49390,49433,49434,49435,49437,49438,49439,49440,49445,49446,49447,
       49448,49449,49451,49452,49453,49454,49455,49456,49457,49458,49464,49466,49467,49650,49651,49652,49653,49655,49656,49658,49660,
       49661,49662,49663,49664,49736,49737,49742,49743,49748,49749,49750,49751,49752,49758,49760,49761,49762,50343,50507,50546,50554,
       50595,50599,50600,50601,50602,50603,50605,50606,50607,50613,50615,50650,50655,50658,50662,50663,50664,50666,50668,50670,50671,
       50673,50674,50676,50677,50678};
    
    
    std::map<int, std::string> triggerNameMap2 = {
        {0, "Clock"},
        {1, "ZDC_South"},
        {2, "ZDC_North"},
        {3, "ZDC_Coincidence"},
        {4, "HCAL_Singles"},
        {5, "HCAL_Coincidence"},
        {8, "MBD_S_geq_1"},
        {9, "MBD_N_geq_1"},
        {10, "MBD_NandS_geq_1"},
        {11, "MBD_NandS_geq_2"},
        {12, "MBD_NandS_geq_1_vtx_lessThen_10_cm"},
        {13, "MBD_NandS_geq_1_vtx_lessThen_30_cm"},
        {14, "MBD_NandS_geq_1_vtx_lessThen_60_cm"},
        {15, "HCAL_Singles_plus_MBD_NS_geq_1"},
        {16, "Jet_6_GeV_plus_MBD_NS_geq_1"},
        {17, "Jet_8_GeV_plus_MBD_NS_geq_1"},
        {18, "Jet_10_GeV_plus_MBD_NS_geq_1"},
        {19, "Jet_12_GeV_plus_MBD_NS_geq_1"},
        {20, "Jet_6_GeV"},
        {21, "Jet_8_GeV"},
        {22, "Jet_10_GeV"},
        {23, "Jet_12_GeV"},
        {24, "Photon_2_GeV_plus_MBD_NS_geq_1"},
        {25, "Photon_3_GeV_plus_MBD_NS_geq_1"},
        {26, "Photon_4_GeV_plus_MBD_NS_geq_1"},
        {27, "Photon_5_GeV_plus_MBD_NS_geq_1"},
        {28, "Photon_2_GeV"},
        {29, "Photon_3_GeV"},
        {30, "Photon_4_GeV"},
        {31, "Photon_5_GeV"}
    };
    
    
    // Determine which map to use
    if (std::find(runNumbersForMap1.begin(), runNumbersForMap1.end(), runNumber) != runNumbersForMap1.end()) {
        std::cout << "Using triggerNameMap1 for run number " << runNumber << std::endl;
        // Use triggerNameMap1
        for (const auto& pair : triggerNameMap1) {
            triggerNameToIndexMap[pair.second] = pair.first;
        }
    } else if (std::find(runNumbersForMap2.begin(), runNumbersForMap2.end(), runNumber) != runNumbersForMap2.end()) {
        std::cout << "Using triggerNameMap2 for run number " << runNumber << std::endl;
        // Use triggerNameMap2
        for (const auto& pair : triggerNameMap2) {
            triggerNameToIndexMap[pair.second] = pair.first;
        }
    } else {
        std::cerr << "Run number " << runNumber << " not found in any trigger map." << std::endl;
    }
}

std::string get_trigger_name_from_histogram(const std::string& histName, const std::map<std::string, int>& triggerNameToIndexMap) {
    // Iterate over all trigger names
    for (const auto& pair : triggerNameToIndexMap) {
        const std::string& triggerName = pair.first;
        // Check if histName ends with triggerName
        if (histName.length() >= triggerName.length() &&
            histName.compare(histName.length() - triggerName.length(), triggerName.length(), triggerName) == 0) {
            return triggerName;
        }
    }
    // If no trigger name found, return empty string
    return "";
}

int scale_histogram(TH1* hist, int scaledown, const std::string& triggerName) {
    if (!hist) {
        std::cerr << ANSI_RED << "Error: Histogram is null for trigger " << triggerName << ANSI_RESET << std::endl;
        return -1; // Error
    }

     // Print the histogram entries, integral, sum of weights
     std::cout << "Histogram \"" << hist->GetName() << "\" has GetEntries() = " << hist->GetEntries()
               << ", Integral() = " << hist->Integral()
               << ", GetSumOfWeights() = " << hist->GetSumOfWeights() << "." << std::endl;

     // Check if histogram has zero content
     if (hist->Integral() == 0) {
         std::cout << ANSI_CYAN << "Skipping histogram \"" << hist->GetName()
                   << "\" for trigger " << triggerName
                   << " because it has zero content (integral is zero)." << ANSI_RESET << std::endl;
         return 1; // Skipped zero content
     }
    
    // Apply the scaling if the scaledown factor is valid
    if (scaledown > 0) {
        hist->Scale(scaledown + 1.0);  // Scale by (scaledown + 1)
        std::cout << ANSI_YELLOW << "Scaled histogram \"" << hist->GetName()
                  << "\" for trigger " << triggerName
                  << " by factor " << scaledown + 1
                  << " (scale-down factor: " << scaledown << ")" << ANSI_RESET << std::endl;
        return 0; // Success
    } else if (scaledown == 0) {
        std::cout << ANSI_CYAN << "No scaling applied to histogram \""
                  << hist->GetName() << "\" for trigger "
                  << triggerName << " (scale-down factor is 0)" << ANSI_RESET << std::endl;
        return 2; // No scaling needed
    } else {
        std::cout << ANSI_RED << "Skipping histogram \"" << hist->GetName()
                  << "\" for trigger " << triggerName
                  << " because the trigger was off (scale-down factor = -1)" << ANSI_RESET << std::endl;
        return 3; // Skipped due to trigger off
    }
}

void processHistogramsInDirectory(TDirectory* dir, const std::map<std::string, int>& triggerNameToIndexMap, const int scaledowns[]) {
    dir->cd();
    TList* keyList = dir->GetListOfKeys();
    if (!keyList) return;

    // Make a copy of the key list to avoid modification during iteration
    std::vector<TKey*> keys;
    TIter next(keyList);
    TKey* key;
    while ((key = (TKey*)next())) {
        keys.push_back((TKey*)key->Clone());
    }

    // Counters
    int histogramsProcessed = 0;
    int histogramsSkippedZeroEntries = 0;
    int histogramsSkippedTriggerOff = 0;
    int histogramsSkippedNoTriggerName = 0;
    int histogramsSkippedNotTH1 = 0;
    int histogramsNoScalingNeeded = 0;

    std::cout << "Processing directory: " << dir->GetPath() << std::endl;

    // Now process each key from the copied list
    for (auto key : keys) {
        dir->cd();
        TObject* obj = key->ReadObj();
        if (!obj) {
            delete key; // Clean up
            continue;
        }

        if (obj->InheritsFrom("TH1")) {
            TH1* hist = (TH1*)obj;
            std::string histName = hist->GetName();

            // Extract trigger name from histogram name
            std::string triggerName = get_trigger_name_from_histogram(histName, triggerNameToIndexMap);

            // Get trigger index
            auto it = triggerNameToIndexMap.find(triggerName);
            if (it == triggerNameToIndexMap.end()) {
                std::cerr << "Trigger name \"" << triggerName << "\" not found in triggerNameToIndexMap" << std::endl;
                histogramsSkippedNoTriggerName++;
                delete hist; // Clean up
                delete key;
                continue;
            }
            int triggerIndex = it->second;

            // Get scaledown factor
            int scaledown = scaledowns[triggerIndex];

            // Scale the histogram
            int scaleResult = scale_histogram(hist, scaledown, triggerName);

            if (scaleResult == 0) {
                histogramsProcessed++;
            } else if (scaleResult == 1) {
                histogramsSkippedZeroEntries++;
            } else if (scaleResult == 3) {
                histogramsSkippedTriggerOff++;
            } else if (scaleResult == 2) {
                histogramsNoScalingNeeded++;
            }

            dir->cd();
            hist->Write(hist->GetName(), TObject::kOverwrite);
            // Delete the histogram object (since we own it)
            delete obj;
            // Delete the key
            delete key;

        } else if (obj->InheritsFrom("TDirectory")) {
            TDirectory* subDir = (TDirectory*)obj;

            // Print statement for entering subdirectory
            std::cout << "Entering subdirectory: " << subDir->GetName() << std::endl;

            processHistogramsInDirectory(subDir, triggerNameToIndexMap, scaledowns);
            dir->cd();
            
            // Do not delete the directory object here; ROOT manages it
            // Only delete the key
            delete key;
        } else {
            histogramsSkippedNotTH1++;
            delete obj;
            delete key;
        }
    }

    // After processing all histograms in this directory, print statistics
    std::cout << "Directory " << dir->GetPath() << " processing summary:" << std::endl;
    std::cout << "  Histograms processed (scaled): " << histogramsProcessed << std::endl;
    std::cout << "  Histograms with no scaling needed: " << histogramsNoScalingNeeded << std::endl;
    std::cout << "  Histograms skipped (zero entries): " << histogramsSkippedZeroEntries << std::endl;
    std::cout << "  Histograms skipped (trigger off): " << histogramsSkippedTriggerOff << std::endl;
    std::cout << "  Histograms skipped (no trigger name): " << histogramsSkippedNoTriggerName << std::endl;
    std::cout << "  Objects skipped (not TH1): " << histogramsSkippedNotTH1 << std::endl;
}


// Function to process QA histograms
void processQAHistograms(const std::string& outputFileName, int runNumber) {
    std::cout << "Starting processing of QA histograms for run " << runNumber << std::endl;

    // Retrieve scale-down factors for the current run
    int scaledowns[64] = {0};
    get_scaledowns(runNumber, scaledowns);

    // Get the trigger name to index map for this run
    std::map<std::string, int> triggerNameToIndexMap;
    get_trigger_name_to_index_map(runNumber, triggerNameToIndexMap);

    // Open the merged ROOT file
    TFile* file = TFile::Open(outputFileName.c_str(), "UPDATE");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open merged ROOT file: " << outputFileName << std::endl;
        return;
    }

    // Process the 'QA' directory
    TDirectory* qaDir = (TDirectory*)file->Get("QA");
    if (!qaDir) {
        std::cerr << "Error: QA directory not found in the ROOT file: " << outputFileName << std::endl;
        file->Close();
        return;
    }

    std::cout << "Processing 'QA' directory" << std::endl;
    processHistogramsInDirectory(qaDir, triggerNameToIndexMap, scaledowns);

    // Similarly process the 'PhotonAnalysis' directory
    TDirectory* photonDir = (TDirectory*)file->Get("PhotonAnalysis");
    if (photonDir) {
        std::cout << "Processing 'PhotonAnalysis' directory" << std::endl;
        processHistogramsInDirectory(photonDir, triggerNameToIndexMap, scaledowns);
    } else {
        std::cerr << "PhotonAnalysis directory not found in the ROOT file: " << outputFileName << std::endl;
    }

    file->Close();
    std::cout << "Successfully processed and scaled histograms for run " << runNumber << std::endl;
}

// Function to merge run files
void mergeRunFiles(const std::string& runNumber, const std::string& baseDir, const std::string& outputDir) {
    std::string runPath = baseDir + runNumber + "/";
    std::string outputFileName = outputDir + runNumber + "_HistOutput.root";

    struct stat buffer;
    if (stat(outputFileName.c_str(), &buffer) == 0) {
        std::cout << "Output file already exists for run " << runNumber << ": " << outputFileName
                  << ". Skipping merge for this run." << std::endl;
        return;
    }

    std::cout << "Processing run number: " << runNumber << std::endl;

    // Open the directory containing the ROOT files
    DIR* dir = opendir(runPath.c_str());
    if (!dir) {
        std::cerr << "Error: Cannot open directory " << runPath << std::endl;
        return;
    }

    std::vector<std::string> rootFiles;
    struct dirent* entry;
    while ((entry = readdir(dir))) {
        std::string fileName = entry->d_name;
        if (fileName.find(".root") != std::string::npos) {
            std::string filePath = runPath + fileName;
            rootFiles.push_back(filePath);
        }
    }
    closedir(dir);

    if (rootFiles.empty()) {
        std::cerr << "No ROOT files found in: " << runPath << std::endl;
        return;
    }
    // Validate each ROOT file before merging
    std::vector<std::string> validRootFiles;
    std::string problematicFileList = "problematicHaddRuns.txt";
    std::ofstream problematicFileStream;

    // Open the problematic file list in append mode
    problematicFileStream.open(problematicFileList, std::ios::app);
    if (!problematicFileStream.is_open()) {
        std::cerr << "Error: Cannot open or create " << problematicFileList << std::endl;
    }

    for (const auto& file : rootFiles) {
        TFile* f = TFile::Open(file.c_str(), "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "Problematic file detected: " << file << std::endl;
            if (problematicFileStream.is_open()) {
                problematicFileStream << runNumber << " " << file << std::endl;
            }
            delete f;
            continue; // Skip this file
        } else {
            validRootFiles.push_back(file);
        }
        delete f;
    }

    if (problematicFileStream.is_open()) {
        problematicFileStream.close();
    }

    if (validRootFiles.empty()) {
        std::cerr << "No valid ROOT files to merge for run " << runNumber << std::endl;
        return;
    }

    // Construct the hadd command to merge all segment files into one output file
    std::string haddCommand = "hadd -f -T " + outputFileName;
    for (const auto& file : rootFiles) {
        haddCommand += " " + file;
    }

    // Execute the hadd command to merge the files
    std::cout << "Merging files for run: " << runNumber << std::endl;
    std::cout << "Executing command: " << haddCommand << std::endl;

    int haddResult = gSystem->Exec(haddCommand.c_str());

    if (haddResult != 0) {
        std::cerr << "Error: hadd failed with exit code " << haddResult << " for run " << runNumber << std::endl;
        return;
    }

    std::cout << "Successfully merged files into " << outputFileName << std::endl;

    // After merging, process and scale histograms in the QA directory
    processQAHistograms(outputFileName, std::stoi(runNumber));
}

// Main function to process a single run number
void mergeSegmentFilesForRuns(Int_t runNumber) {
    std::cout << "Starting mergeSegmentFilesForRuns for run number " << runNumber << std::endl;

    std::string runNumberStr = std::to_string(runNumber);
    std::string baseDir = "/sphenix/tg/tg01/bulk/jbennett/DirectPhotons/output/";
    std::string outputDir = "/sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/output/";

    // Merge the run files and process histograms
    mergeRunFiles(runNumberStr, baseDir, outputDir);
}

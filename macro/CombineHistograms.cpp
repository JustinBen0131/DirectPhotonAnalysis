#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TDirectory.h>
#include <TString.h>
#include <TSystem.h>
#include <iostream>
#include <vector>
#include <string>
#include <dirent.h>
#include <iomanip>

// ANSI escape codes for formatting
#define RESET "\033[0m"
#define BOLD "\033[1m"
#define RED "\033[31m"
#define GREEN "\033[32m"

// Set constant column widths for better alignment
const int colWidthName = 50;
const int colWidthCount = 20;

void get_scaledowns(int runnumber, int scaledowns[]) {
  TSQLServer *db = TSQLServer::Connect("pgsql://sphnxdaqdbreplica:5432/daq","phnxro","");
 
  if (db) {
      printf("Server info: %s\n", db->ServerInfo());
  }
  else {
      printf("bad\n");
  }
  TSQLRow *row;
  TSQLResult *res;
  TString cmd = "";
  char sql[1000];
  
  for (int is = 0; is < 64; is++) {
      sprintf(sql, "select scaledown%02d from gl1_scaledown where runnumber = %d;", is, runnumber);
      printf("%s \n" , sql);
 
      res = db->Query(sql);
      
      int nrows = res->GetRowCount();

      int nfields = res->GetFieldCount();
      for (int i = 0; i < nrows; i++) {
          row = res->Next();
          for (int j = 0; j < nfields; j++) {
              scaledowns[is] = stoi(row->GetField(j));
          }
          delete row;
      }
      delete res;
  }
  delete db;
}

void CombineHistograms() {
    // Input and output directories
    TString inputDir = "/sphenix/tg/tg01/bulk/jbennett/DirectPhotons/outputHists";
    TString outputDir = "/sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/output";
    
    std::cout << "Enters input directory: " << inputDir << std::endl;
    
    // Open the input directory
    DIR* dir = opendir(inputDir.Data());
    if (!dir) {
        std::cerr << "Cannot open input directory: " << inputDir << std::endl;
        return;
    }
    
    struct dirent* entry;
    std::vector<TString> runNumbers;
    
    // Read all run numbers (folders) in the directory
    while ((entry = readdir(dir)) != nullptr) {
        if (entry->d_type == DT_DIR) {
            TString runNumber(entry->d_name);
            if (runNumber.IsDigit()) { // Ensures the folder name is a valid run number
                runNumbers.push_back(runNumber);
            }
        }
    }
    closedir(dir);
    
    std::cout << BOLD << "Run numbers to loop over include: " << RESET;
    for (const auto& runNumber : runNumbers) {
        std::cout << runNumber << " ";
    }
    std::cout << "\n--------------------------------------------------------------------------" << std::endl;
    
    std::vector<TString> combinedRunFiles; // Store the combined files from each run
    
    // Process each run number
    for (const TString& runNumber : runNumbers) {
        TString runDir = inputDir + "/" + runNumber;
        std::cout << BOLD << "Run number = " << runNumber << RESET << std::endl;
        
        // Get list of ROOT files (segments) in the run directory
        DIR* runDirHandle = opendir(runDir.Data());
        if (!runDirHandle) {
            std::cerr << "Cannot open run number directory: " << runDir << std::endl;
            continue;
        }
        
        struct dirent* fileEntry;
        std::vector<TString> rootFiles;
        while ((fileEntry = readdir(runDirHandle)) != nullptr) {
            TString fileName(fileEntry->d_name);
            if (fileName.EndsWith(".root")) {
                rootFiles.push_back(runDir + "/" + fileName);
            }
        }
        closedir(runDirHandle);
        
        std::cout << "Number of segment ROOT files to add = " << rootFiles.size() << std::endl;
        
        // Create output ROOT file for this run number
        TString outputFileName = outputDir + "/combined_" + runNumber + ".root";
        TFile* outputFile = new TFile(outputFileName, "RECREATE");
        
        combinedRunFiles.push_back(outputFileName);  // Add this file to the list of run outputs
        
        // Open the first ROOT file and copy directory structure and histograms
        TFile* firstFile = TFile::Open(rootFiles[0]);
        if (!firstFile || firstFile->IsZombie()) {
            std::cerr << "Cannot open first file: " << rootFiles[0] << std::endl;
            delete outputFile;
            continue;
        }
        
        std::cout << "Entering root file: " << rootFiles[0] << std::endl;
        
        // Clone directory structure and copy histograms from the first file
        TDirectory* invMassDir = (TDirectory*)firstFile->Get("InvariantMassDistributions");
        TDirectory* qaDir = (TDirectory*)firstFile->Get("QA");
        
        if (!invMassDir || !qaDir) {
            std::cerr << "Missing required directories in first file!" << std::endl;
            firstFile->Close();
            delete outputFile;
            continue;
        }
        
        outputFile->cd();
        TDirectory* invMassOutDir = outputFile->mkdir("InvariantMassDistributions");
        TDirectory* qaOutDir = outputFile->mkdir("QA");
        
        // Copy histograms to the new output file
        invMassDir->cd();
        invMassOutDir->cd();
        TIter next(invMassDir->GetListOfKeys());
        TKey* key;
        
        std::cout << BOLD << "Found histograms in InvariantMassDistributions subdirectory:" << RESET << std::endl;
        std::cout << std::setw(colWidthName) << BOLD << "Hist names" << std::setw(colWidthCount) << "Event count" << RESET << std::endl;
        
        while ((key = (TKey*)next())) {
            TString histName = key->GetName();
            TH1* histIn = (TH1*)key->ReadObj();
            std::cout << std::setw(colWidthName) << histName << std::setw(colWidthCount) << histIn->GetEntries() << std::endl;
            
            TH1* histOut = (TH1*)invMassOutDir->Get(histName);
            if (histOut) {
                histOut->Add(histIn);
            } else {
                invMassOutDir->cd();
                histIn->Clone(histName);
            }
        }
        
        // Process the rest of the files and update the output file
        for (size_t i = 1; i < rootFiles.size(); ++i) {
            TFile* inFile = TFile::Open(rootFiles[i]);
            if (!inFile || inFile->IsZombie()) {
                std::cerr << "Cannot open file: " << rootFiles[i] << std::endl;
                continue;
            }
            
            std::cout << "Entering root file: " << rootFiles[i] << std::endl;
            
            // Invariant Mass Directory
            invMassDir = (TDirectory*)inFile->Get("InvariantMassDistributions");
            invMassOutDir->cd();
            if (invMassDir) {
                invMassDir->cd();
                TIter next(invMassDir->GetListOfKeys());
                
                std::cout << BOLD << "Found histograms in InvariantMassDistributions subdirectory:" << RESET << std::endl;
                std::cout << std::setw(colWidthName) << BOLD << "Hist names" << std::setw(colWidthCount) << "Event count" << RESET << std::endl;
                
                while ((key = (TKey*)next())) {
                    TString histName = key->GetName();
                    TH1* histIn = (TH1*)key->ReadObj();
                    std::cout << std::setw(colWidthName) << histName << std::setw(colWidthCount) << histIn->GetEntries() << std::endl;
                    
                    TH1* histOut = (TH1*)invMassOutDir->Get(histName);
                    if (histOut) {
                        histOut->Add(histIn);
                        std::cout << RED << "Added to output: " << std::setw(colWidthName) << histName << std::setw(colWidthCount) << histOut->GetEntries() << RESET << std::endl;
                    } else {
                        invMassOutDir->cd();
                        histIn->Clone(histName);
                        std::cout << RED << "Added to output (new): " << std::setw(colWidthName) << histName << std::setw(colWidthCount) << histIn->GetEntries() << RESET << std::endl;
                    }
                }
            }
            
            // QA Directory
            qaDir = (TDirectory*)inFile->Get("QA");
            qaOutDir->cd();
            if (qaDir) {
                qaDir->cd();
                TIter nextQa(qaDir->GetListOfKeys());
                TKey* qaKey;
                
                std::cout << BOLD << "Found histograms in QA subdirectory:" << RESET << std::endl;
                std::cout << std::setw(colWidthName) << BOLD << "Hist names" << std::setw(colWidthCount) << "Event count" << RESET << std::endl;
                
                while ((qaKey = (TKey*)nextQa())) {
                    TString qaHistName = qaKey->GetName();
                    TH1* qaHistIn = (TH1*)qaKey->ReadObj();
                    
                    TH1* qaHistOut = (TH1*)qaOutDir->Get(qaHistName);
                    if (qaHistOut) {
                        qaHistOut->Add(qaHistIn);
                        std::cout << RED << "Added to output: " << std::setw(colWidthName) << qaHistName << std::setw(colWidthCount) << qaHistOut->GetEntries() << RESET << std::endl;
                    } else {
                        qaOutDir->cd();
                        qaHistIn->Clone(qaHistName);
                        std::cout << RED << "Added to output (new): " << std::setw(colWidthName) << qaHistName << std::setw(colWidthCount) << qaHistIn->GetEntries() << RESET << std::endl;
                    }
                }
            }
            
            inFile->Close();
            std::cout << GREEN << "All histograms from " << rootFiles[i] << " added to output ROOT file." << RESET << std::endl;
        }
        
        // Print final tally for the output root file
        std::cout << RED << BOLD << "Final output root file has histograms:" << RESET << std::endl;
        invMassOutDir->cd();
        TIter nextOut(invMassOutDir->GetListOfKeys());
        TKey* keyOut;
        std::cout << std::setw(colWidthName) << BOLD << "Hist names" << std::setw(colWidthCount) << "Event count" << RESET << std::endl;
        while ((keyOut = (TKey*)nextOut())) {
            TH1* histOut = (TH1*)keyOut->ReadObj();
            std::cout << std::setw(colWidthName) << histOut->GetName() << std::setw(colWidthCount) << histOut->GetEntries() << std::endl;
        }
        
        outputFile->Write();
        outputFile->Close();
        std::cout << RED << "Output root file saved to: " << outputFileName << RESET << std::endl;
        std::cout << "-----------------------------------------------------------------------------------------------" << std::endl;
    }
    
    std::cout << "All run numbers processed!" << std::endl;
    
    // Now, process all the combined run files to generate AllRuns_histograms.root
    TString finalOutputFileName = outputDir + "/AllRuns_histograms.root";
    TFile* finalOutputFile = new TFile(finalOutputFileName, "RECREATE");
    
    TDirectory* finalInvMassDir = finalOutputFile->mkdir("InvariantMassDistributions");
    TDirectory* finalQADir = finalOutputFile->mkdir("QA");
    
    // Iterate through each combined run file and aggregate histograms
    for (const TString& combinedFile : combinedRunFiles) {
        TFile* runFile = TFile::Open(combinedFile);
        if (!runFile || runFile->IsZombie()) {
            std::cerr << "Cannot open combined file: " << combinedFile << std::endl;
            continue;
        }
        
        std::cout << GREEN << "Processing combined file: " << combinedFile << RESET << std::endl;
        
        TDirectory* invMassDir = (TDirectory*)runFile->Get("InvariantMassDistributions");
        TDirectory* qaDir = (TDirectory*)runFile->Get("QA");
        
        int nEvents = 1;  
        
        TH1* eventCountHist = (TH1*)runFile->Get("EventCount");
        if (eventCountHist) {
            nEvents = eventCountHist->GetEntries();
            std::cout << "Number of events in the file: " << nEvents << std::endl;
        }
        
        if (invMassDir) {
            finalInvMassDir->cd();
            TIter nextHist(invMassDir->GetListOfKeys());
            TKey* key;
            while ((key = (TKey*)nextHist())) {
                TString histName = key->GetName();
                TH1* histIn = (TH1*)key->ReadObj();
                TH1* histOut = (TH1*)finalInvMassDir->Get(histName);
                if (histOut) {
                    histOut->Add(histIn);
                } else {
                    finalInvMassDir->cd();
                    histIn->Clone(histName);
                }
            }
        }
        
        std::vector<TString> histogramsToNormalize = {
            "hTotalCaloEEMCal",
            "hTotalCaloEOHCal",
            "hTotalCaloEIHCal",
            "hClusterChi2",
            "h_ohcalChi2",
            "h_ihcalChi2",
            "hClusterECore",
            "hClusterPt",
            "hVtxZ"
        };
        
        if (qaDir) {
            finalQADir->cd();
            TIter nextQaHist(qaDir->GetListOfKeys());
            TKey* qaKey;
            while ((qaKey = (TKey*)nextQaHist())) {
                TString histName = qaKey->GetName();
                TH1* qaHistIn = (TH1*)qaKey->ReadObj();
                
                // Check if the histogram is in the list of histograms to normalize
                bool shouldNormalize = std::find(histogramsToNormalize.begin(), histogramsToNormalize.end(), histName) != histogramsToNormalize.end();
                
                // Normalize the QA histogram by 1/nEvents only if it's in the list
                if (shouldNormalize && qaHistIn->GetEntries() > 0) {
                    qaHistIn->Scale(1.0 / nEvents);
                    std::cout << "Normalized histogram: " << histName << " by 1/nEvents" << std::endl;
                }
                
                TH1* qaHistOut = (TH1*)finalQADir->Get(histName);
                if (qaHistOut) {
                    qaHistOut->Add(qaHistIn);
                    std::cout << RED << "Added to output: " << std::setw(colWidthName) << histName << std::setw(colWidthCount) << qaHistOut->GetEntries() << RESET << std::endl;
                } else {
                    finalQADir->cd();
                    qaHistIn->Clone(histName);
                    std::cout << RED << "Added to output (new): " << std::setw(colWidthName) << histName << std::setw(colWidthCount) << qaHistIn->GetEntries() << RESET << std::endl;
                }
            }
        }
        
        runFile->Close();
        std::cout << GREEN << "All histograms from " << combinedFile << " added to final output ROOT file." << RESET << std::endl;
    }
    
    // Save and close the final output file
    finalOutputFile->Write();
    finalOutputFile->Close();
    
    std::cout << RED << "Final output root file for all runs saved to: " << finalOutputFileName << RESET << std::endl;
}

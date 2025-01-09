//////////////////////////////////////////////////////////////////////////
//
//  SplitCaloAna.C
//
//  PURPOSE:
//    (1) Open a large ROOT file containing TTree "slimtree"
//    (2) Split it into multiple "segment" files in 200k-event chunks
//    (3) Write a ".list" file of the split segment paths
//    (4) Ignore all other objects (like histograms) in the parent file
//
//  MAJOR SIMPLIFICATION:
//    Use TTree::CopyTree(..., nentries, firstentry) to copy partial
//    ranges of the TTree, letting ROOT handle the branch setup
//    and reading/filling.
//
//  USAGE:
//    root -l -b -q SplitCaloAna.C
//
//////////////////////////////////////////////////////////////////////////

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <string>

void SplitCaloAna()
{
  //////////////////////////////////////////////////////////////////////////
  // 1) OPEN INPUT FILE
  //////////////////////////////////////////////////////////////////////////
  const TString inFileName =
    "/sphenix/user/patsfan753/tutorials/tutorials/"
    "CaloDataAnaRun24pp/SimulationPhoton10_Jan7_25/caloana.root";

  std::cout << "[INFO] Opening input file: " << inFileName << std::endl;
  TFile* inFile = TFile::Open(inFileName, "READ");
  if (!inFile || inFile->IsZombie())
  {
    std::cerr << "[ERROR] Could not open input file: " << inFileName << "\n";
    return;
  }

  // Grab TTree "slimtree"
  TTree* inTree = dynamic_cast<TTree*>( inFile->Get("slimtree") );
  if (!inTree)
  {
    std::cerr << "[ERROR] 'slimtree' not found in the input file.\n";
    inFile->Close();
    return;
  }
  const Long64_t nEntries = inTree->GetEntries();
  std::cout << "[INFO] 'slimtree' has " << nEntries << " entries.\n";


  //////////////////////////////////////////////////////////////////////////
  // 2) DETERMINE CHUNK SIZES
  //////////////////////////////////////////////////////////////////////////
  const Long64_t chunkSize = 200000; // 200k events per segment
  Long64_t nSegments = nEntries / chunkSize;
  if ((nEntries % chunkSize) != 0)
    nSegments++;

  std::cout << "[INFO] Will create " << nSegments
            << " segment(s). Each has " << chunkSize
            << " events, except possibly the last leftover.\n";


  //////////////////////////////////////////////////////////////////////////
  // 3) PREPARE OUTPUT DIRECTORY
  //////////////////////////////////////////////////////////////////////////
  const TString outDir =
    "/sphenix/u/patsfan753/scratch/tutorials/tutorials/"
    "CaloDataAnaRun24pp/SimulationPhoton10_Jan7_25/caloanaSegments";

  // Ensure it exists
  if (gSystem->mkdir(outDir, kTRUE) != 0 && gSystem->AccessPathName(outDir, kFileExists))
  {
    std::cerr << "[ERROR] Cannot create output directory: " << outDir << "\n";
    inFile->Close();
    return;
  }


  //////////////////////////////////////////////////////////////////////////
  // 4) OPEN .list FILE FOR WRITING ALL SEGMENT NAMES
  //////////////////////////////////////////////////////////////////////////
  const TString listFilePath =
    "/sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/"
    "caloAnaSimListFile.list";

  std::ofstream listFile(listFilePath.Data());
  if (!listFile.is_open())
  {
    std::cerr << "[WARNING] Could not open .list file:\n  "
              << listFilePath << "\n(No segment paths will be written.)\n";
  }


  //////////////////////////////////////////////////////////////////////////
  // 5) SPLIT LOOP OVER SEGMENTS
  //////////////////////////////////////////////////////////////////////////
  Long64_t currentStart = 0;
  for (Long64_t seg = 1; seg <= nSegments; seg++)
  {
    // Compute event range for this segment
    Long64_t startEvt = currentStart;
    Long64_t endEvt   = startEvt + chunkSize;
    if (endEvt > nEntries) endEvt = nEntries;
    Long64_t numThisChunk = endEvt - startEvt;

    std::cout << "\n[INFO] Segment #" << seg << " of " << nSegments
              << ": events [" << startEvt << ", " << (endEvt - 1)
              << "], " << numThisChunk << " total.\n";

    // Construct output filename
    TString outFileName = Form("%s/caloAnaSegments_%lld.root", outDir.Data(), seg);

    // Create the output file
    TFile* outFile = TFile::Open(outFileName, "RECREATE");
    if (!outFile || outFile->IsZombie())
    {
      std::cerr << "[ERROR] Could not create: " << outFileName << "\n";
      if (outFile) { outFile->Close(); delete outFile; }
      currentStart = endEvt;
      continue;
    }

    //---------------------------------------------------------
    // Copy the sub-range of "slimtree" only
    //---------------------------------------------------------
    outFile->cd();
    TTree* outTree = inTree->CopyTree("", "", numThisChunk, startEvt);
    if (!outTree)
    {
      std::cerr << "[ERROR] TTree->CopyTree(...) returned null!\n";
    }
    else
    {
      outTree->SetName("slimtree");
      outTree->Write("slimtree");
      std::cout << "[INFO] segment #" << seg << " => " << outFileName
                << " with " << outTree->GetEntries() << " event(s).\n";
    }

    // Finalize
    outFile->Write();
    outFile->Close();
    delete outFile;
    outFile = nullptr;

    //---------------------------------------------------------
    // Optional post-check (verifying TTree entry count)
    //---------------------------------------------------------
    {
      TFile* checkFile = TFile::Open(outFileName, "READ");
      if (!checkFile || checkFile->IsZombie())
      {
        std::cerr << "[WARNING] Could not re-open segment for verify: "
                  << outFileName << "\n";
        if (checkFile) delete checkFile;
      }
      else
      {
        TTree* verifyTree = dynamic_cast<TTree*>( checkFile->Get("slimtree") );
        if (!verifyTree)
        {
          std::cerr << "[ERROR] 'slimtree' missing in " << outFileName << "!\n";
        }
        else
        {
          std::cout << "[INFO] Verified segment has "
                    << verifyTree->GetEntries() << " events.\n";
        }
        checkFile->Close();
        delete checkFile;
      }
    }

    //---------------------------------------------------------
    // Write segment path to the .list file
    //---------------------------------------------------------
    if (listFile.is_open()) {
      listFile << outFileName << "\n";
    }

    // Move to next chunk
    currentStart = endEvt;
  } // end for (seg)


  //////////////////////////////////////////////////////////////////////////
  // 6) CLEANUP
  //////////////////////////////////////////////////////////////////////////
  inFile->Close();
  delete inFile;

  if (listFile.is_open())
  {
    listFile.close();
    std::cout << "\n[INFO] Wrote segment list to:\n  " << listFilePath << "\n";
  }
}

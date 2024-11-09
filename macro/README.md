mergeSegmentFilesForRuns.C adds segments from analysis output for run by run root files, querying scaledown from sphenix db and scaling necessary histograms

runCondor_hAddRunByRunSegments.sh -- running ./runCondor_hAddRunByRunSegments.sh submits condor job per run number to scale and add run by run output in //

mergeSegmentFilesForRuns.sh is the executable for ^

mergeRunFilesToFinalOutput.C reads in run by run root files and hadds for final root file

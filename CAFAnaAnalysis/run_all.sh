cafe -bq Scripts/Selection.cxx
cafe -bq Scripts/SelectionInteBreakdown.C
cafe -bq Scripts/SelectionTopologyBreakdown.C
cafe -bq Scripts/SelectionCutPlots.C
root -b -q Scripts/SerialPlotGenerator.cpp
cafe -bq Scripts/SelectionRunData.C
cafe -bq Scripts/SelectionEfficiency.C
cafe -bq Scripts/SelectionMigrationMatrix.C
cafe -bq Scripts/RunAllSystematics.C
cafe -bq Scripts/SelectionNTargetSystematics.C
cafe -bq Scripts/SelectionPOTSystematics.C
cafe -bq Scripts/SelectionReinteractionSystematics.C
cafe -bq Scripts/SelectionDetectorSystematics.C
cafe -bq Scripts/SelectionMCStatSystematics.C
root -b -q Scripts/StatSystematics.cpp
root -b -q Scripts/TotalCovMatrices.cpp
root -b -q Scripts/Unfold.cpp
root -b -q Scripts/WienerSVDOverlay.cpp
cafe -bq Scripts/SelectionFakeData.C
root -b -q Scripts/UnfoldFakeData.cpp
cafe -bq Scripts/SelectionData.C

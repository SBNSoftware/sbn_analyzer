{

    vector<TString> WhichSample; vector<TString> WhichName;

    //----------------------------------------//

    // Normal samples
    WhichSample.push_back("bnb.ub.num.nuwro_19_02_1.flat"); WhichName.push_back("NuWro");
    WhichSample.push_back("bnb.ub.num.neut_5_4_0_1.flat"); WhichName.push_back("NEUT");
    WhichSample.push_back("NEUT.flat"); WhichName.push_back("NEUT");
    WhichSample.push_back("14_1000180400_CC_v3_4_0_G18_10a_02_11a.flat"); WhichName.push_back("GENIE_G18");
    WhichSample.push_back("14_1000180400_CC_v3_4_0_AR23_20i_00_000.flat"); WhichName.push_back("GENIE_AR23");
    WhichSample.push_back("GiBUU.flat"); WhichName.push_back("GiBUU");

    // Samples without final state interactions
    WhichSample.push_back("GiBUU_noFSI.flat"); WhichName.push_back("GiBUU_NoFSI");

    // Samples with only MEC events
    WhichSample.push_back("14_1000180400_CCMEC_v3_4_2_AR23_20i_00_000.flat"); WhichName.push_back("GENIE_AR23_MEC");
    WhichSample.push_back("14_1000180400_CCMEC_v3_4_2_G18_02a_00_000.flat"); WhichName.push_back("GENIE_G18_Emp_MEC");
    WhichSample.push_back("14_1000180400_CCMEC_v3_4_2_G18_10a_00_000.flat"); WhichName.push_back("GENIE_G18_Nie_MEC");

    //----------------------------------------//

    gROOT->ProcessLine(".L ./../Utils/Tools.cxx+");
    gROOT->ProcessLine(".L ./../Utils/TwoPTools.cxx+");
    gROOT->ProcessLine(".L ./Scripts/FlatTreeAnalyzer.cxx+");

    for (int i =0;i < (int)(WhichSample.size()); i++) {
        gROOT->ProcessLine("FlatTreeAnalyzer(\""+WhichSample[i]+"\",\""+WhichName[i]+"\").Loop()");
    }
    gROOT->ProcessLine(".q");
};

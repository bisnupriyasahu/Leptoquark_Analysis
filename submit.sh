./rootcom LQMETClass all

#CRAB_DATASET1 

#./all output_DYJETS_HT100200.root /hdfs/store/user/varuns/NTuples/MC/MC2018_RunIIAutumn18_JECv19/DYJets/DYJetsToLL_HT100-200/0000/*.root
#./all output_DYJETS_HT100-200_1.root /hdfs/store/user/varuns/NTuples/MC/MC2018_RunIIAutumn18_JECv19/DYJets/DYJetsToLL_HT100-200/0001/MC_DYJetsToLL_HT10*.root
#./all output_DYJETS_HT200-400_0.root /hdfs/store/user/varuns/NTuples/MC/MC2018_RunIIAutumn18_JECv19/DYJets/DYJetsToLL_HT200-400/0000/*.r
#./all TRY_output_DYJETHT2500toinf_1.root DYJETS_HT200-400_0001 /hdfs/store/user/varuns/NTuples/MC/MC2018_RunIIAutumn18_JECv19/DYJets/DYJetsToLL_HT2500-Inf/0000/MC_DYJetsToLL_HT2500-Inf_21.root
#./all TRY_with8rootfiles.root DYJETS_HT200-400_0001 /nfs_scratch/bsahu/LQ_Analysis/CMSSW_10_2_18/src/MuonAnalyzer/Muon_Analyzer/test/rootfile/wjets_rottfiles/*.root
#./all TRY_with8rootfiles.root DYJETS_HT200-400_0001 /hdfs/store/user/varuns/NTuples/MC/MC2018_RunIIAutumn18_JECv19/DYJets/DYJetsToLL_HT70-100/0000/*1*.root



#CRAB_DATASET_CRAB_DYJetsToLL_HT70-100

#./MakeCondorFiles.sh all output_DYJETS_HT70-100_0000_noplotfill_down.root DYJETS_HT70-100_0000_noplotfill_down /hdfs/store/user/varuns/NTuples/MC/MC2018_RunIIAutumn18_JECv19/DYJets/DYJetsToLL_HT70-100/0000/
#./MakeCondorFiles.sh all output_DYJETS_HT70-100_0000_noplotfill_up_down.root DYJETS_HT70-100_0000_noplotfill_up_down /hdfs/store/user/varuns/NTuples/MC/MC2018_RunIIAutumn18_JECv19/DYJets/DYJetsToLL_HT70-100/0000/
./MakeCondorFiles.sh all output_DYJETS_HT70-100_0000_noplotfill.root DYJETS_HT70-100_0000_noplotfill /hdfs/store/user/varuns/NTuples/MC/MC2018_RunIIAutumn18_JECv19/DYJets/DYJetsToLL_HT70-100/0000/

#./MakeCondorFiles.sh all CHK_output_DYJETS_HT70-100_0000_noplotfill.root CHK_DYJETS_HT70-100_0000_noplotfill /hdfs/store/user/varuns/NTuples/MC/MC2018_RunIIAutumn18_JECv19/DYJets/DYJetsToLL_HT70-100/0000/
#./MakeCondorFiles.sh all output_DYJETS_HT70-100_0000.root DYJETS_HT70-100_0000 /hdfs/store/user/varuns/NTuples/MC/MC2018_RunIIAutumn18_JECv19/DYJets/DYJetsToLL_HT70-100/0000/

#./MakeCondorFiles.sh all output_DYJETS_HT70-100_0001.root DYJETS_HT70-100_0001 /hdfs/store/user/varuns/NTuples/MC/MC2018_RunIIAutumn18_JECv19/DYJets/DYJetsToLL_HT70-100/0001/


#CRAB_DATASET_CRAB_DYJetsToLL_HT100-200

#./MakeCondorFiles.sh all output_DYJETS_HT100-200_0000.root DYJETS_HT100-200_0000 /hdfs/store/user/varuns/NTuples/MC/MC2018_RunIIAutumn18_JECv19/DYJets/DYJetsToLL_HT100-200/0000/
#./MakeCondorFiles.sh all output_DYJETS_HT100-200_0001.root DYJETS_HT100-200_0001 /hdfs/store/user/varuns/NTuples/MC/MC2018_RunIIAutumn18_JECv19/DYJets/DYJetsToLL_HT100-200/0001/




#CRAB_DATASET_CRAB_DYJetsToLL_HT200-400

#./MakeCondorFiles.sh all output_DYJETS_HT200-400_0000.root DYJETS_HT200-400_0000 /hdfs/store/user/varuns/NTuples/MC/MC2018_RunIIAutumn18_JECv19/DYJets/DYJetsToLL_HT200-400/0000/
#./MakeCondorFiles.sh all output_DYJETS_HT200-400_0001.root DYJETS_HT200-400_0001_1 /hdfs/store/user/varuns/NTuples/MC/MC2018_RunIIAutumn18_JECv19/DYJets/DYJetsToLL_HT200-400/0001/


#CRAB_DATASET_CRAB_DYJetsToLL_HT400-600
#./MakeCondorFiles.sh all output_DYJETS_HT400-600_0000.root DYJETS_HT400-600_0000 /hdfs/store/user/varuns/NTuples/MC/MC2018_RunIIAutumn18_JECv19/DYJets/DYJetsToLL_HT400-600/0000/
#./MakeCondorFiles.sh all output_DYJETS_HT400-600_0001.root DYJETS_HT400-600_0001 /hdfs/store/user/varuns/NTuples/MC/MC2018_RunIIAutumn18_JECv19/DYJets/DYJetsToLL_HT400-600/0001/



#CRAB_DATASET_CRAB_DYJetsToLL_HT600-800
#./MakeCondorFiles.sh all output_DYJETS_HT600-800_0000.root DYJETS_HT600-800_0000 /hdfs/store/user/varuns/NTuples/MC/MC2018_RunIIAutumn18_JECv19/DYJets/DYJetsToLL_HT600-800/0000/




#CRAB_DATASET_CRAB_DYJetsToLL_HT800-1200
#./MakeCondorFiles.sh all output_DYJETS_HT800-1200_0000.root DYJETS_HT800-1200_0000 /hdfs/store/user/varuns/NTuples/MC/MC2018_RunIIAutumn18_JECv19/DYJets/DYJetsToLL_HT800-1200/0000/
                                                                                                                                                                                 



#CRAB_DATASET_CRAB_DYJetsToLL_HT1200-2500
#./MakeCondorFiles.sh all output_DYJETS_HT1200-2500_0000.root DYJETS_HT1200-2500_0000 /hdfs/store/user/varuns/NTuples/MC/MC2018_RunIIAutumn18_JECv19/DYJets/DYJetsToLL_HT1200-2500/0000/



#CRAB_DATASET_CRAB_DYJetsToLL_HT2500-inf
#./MakeCondorFiles.sh all output_DYJETS_HT2500-Inf_0000.root DYJETS_HT2500-Inf_0000 /hdfs/store/user/varuns/NTuples/MC/MC2018_RunIIAutumn18_JECv19/DYJets/DYJetsToLL_HT2500-Inf/0000/

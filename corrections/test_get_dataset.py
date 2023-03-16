from gbr_trainer import GBRLikelihoodTrainer
gbr = GBRLikelihoodTrainer("params_gbrtrain.txt")
#gbr.get_dataset("~/cernbox/CYGNUS/reco/lime/merged_winter22/reco_runs5X.root","~/cernbox/CYGNUS/reco/lime/merged_winter22/friends/reco_runs5X_Friend.root",0,100)
gbr.get_dataset("~/cernbox/CYGNUS/reco/lime/merged_winter23/merged_feruns_8882_9857.root","~/cernbox/CYGNUS/reco/lime/merged_winter23/merged_feruns_8882_9857_Friend.root",savepanda="test.pkl")

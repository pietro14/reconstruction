# USAGE: python gbr_trainer.py ~/Work/data/cygnus/RECO/lime2021/v2/fe55/reco_run04455_3D.root params_gbrtrain.txt

import ROOT
ROOT.gROOT.SetBatch(True)

import numpy as np
from root_numpy import tree2array

import matplotlib.pyplot as plt
import numpy as np

from sklearn import ensemble
from sklearn.inspection import permutation_importance
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split


class GBRLikelihoodTrainer:
    def set_defaults(self,params):
        self.tree_name = params['tree_name']
        self.target = params['target']
        self.var = params['inputs']
        self.cuts_base = params['selection']
        
    def __init__(self,params):
        self.set_defaults(params)
        training_keys = ['n_estimators','max_depth','min_samples_split','learning_rate']
        self.training_params = {k: params[k] for k in training_keys}

    def name(self):
        return "{args.base_name}_{args.vars_name}_{args.cuts_name}_vgem{args.vgem1}V_ntrees{args.ntrees}".format(args=self)

    def variables(self):
        return self.var.split(":")
    
    def get_dataset(self,rfile):
        tfile = ROOT.TFile.Open(rfile)
        tree = tfile.Get(self.tree_name)

        # so target is always the first variable
        variables = [self.target] + self.var.split(":")
        print("List of variables = ",variables)
        dataset = tree2array(tree,variables,object_selection={self.cuts_base : variables})
        tfile.Close()
        
        conc = np.stack(dataset[0],axis=-1)
        for i in range(1,len(dataset)):
            conc = np.append(conc,np.stack(dataset[i],axis=-1),axis=0)
        X = conc[:,1:]
        y = conc[:,0]
        return X,y

    def train_model(self,X,y):
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=13)
        reg = ensemble.GradientBoostingRegressor(**self.training_params)
        reg.fit(X_train, y_train)

        mse = mean_squared_error(y_test, reg.predict(X_test))
        print("The mean squared error (MSE) on test set: {:.4f}".format(mse))

        self.X_test = X_test
        self.y_test = y_test
        self.reg_ = reg

    def plot_training(self):
        test_score = np.zeros((self.training_params['n_estimators'],), dtype=np.float64)
        for i, y_pred in enumerate(self.reg_.staged_predict(self.X_test)):
            test_score[i] = self.reg_.loss_(self.y_test, y_pred)

        fig = plt.figure(figsize=(6, 6))
        plt.subplot(1, 1, 1)
        plt.title('Deviance')
        plt.plot(np.arange(self.training_params['n_estimators']) + 1, self.reg_.train_score_, 'b-',
                 label='Training Set Deviance')
        plt.plot(np.arange(self.training_params['n_estimators']) + 1, test_score, 'r-',
                 label='Test Set Deviance')
        plt.legend(loc='upper right')
        plt.xlabel('Boosting Iterations')
        plt.ylabel('Deviance')
        plt.savefig('training.png',bbox_inches='tight', pad_inches=0)
      
        # plot variables importance
        feature_importance = self.reg_.feature_importances_
        sorted_idx = np.argsort(feature_importance)
        pos = np.arange(sorted_idx.shape[0]) + .5
        fig = plt.figure(figsize=(12, 6))
        plt.subplot(1, 2, 1)
        plt.barh(pos, feature_importance[sorted_idx], align='center')
        plt.yticks(pos, np.array(self.variables())[sorted_idx])
        plt.title('Feature Importance (MDI)')
        
        result = permutation_importance(self.reg_, self.X_test, self.y_test, n_repeats=10,
                                        random_state=42, n_jobs=2)
        sorted_idx = result.importances_mean.argsort()
        plt.subplot(1, 2, 2)
        plt.boxplot(result.importances[sorted_idx].T,
                    vert=False, labels=np.array(self.variables())[sorted_idx])
        plt.title("Permutation Importance (test set)")
        plt.savefig('variables_importance.png')

        
    def save_model(self,filename):
        import joblib
        joblib.dump(self.reg_,filename)
    
        
if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser(usage='%prog input.root params_gbrtrain.txt [opts] ')
    parser.add_option('-o', '--outname', dest='outname', default='gbrLikelihood.sav', type='string', help='output name of the regression model')
    (options, args) = parser.parse_args()

    recofile = args[0]
    config = open(args[1], "r")
    params = eval(config.read())

    GBR = GBRLikelihoodTrainer(params)
    X,y = GBR.get_dataset(recofile)
    
    GBR.train_model(X,y)
    GBR.plot_training()
    GBR.save_model(options.outname)


# USAGE: python gbr_trainer.py ~/Work/data/cygnus/RECO/lime2021/v2/fe55/reco_run04455_3D.root params_gbrtrain.txt

import ROOT
ROOT.gROOT.SetBatch(True)

import numpy as np
from root_numpy import tree2array,fill_hist

import matplotlib.pyplot as plt
import numpy as np
import joblib

from sklearn import ensemble
from sklearn.inspection import permutation_importance
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split

def getCanvas(name='c'):

    ROOT.gStyle.SetPalette(ROOT.kRainBow)
    ROOT.gStyle.SetNumberContours(51)
    ROOT.gErrorIgnoreLevel = 100
    ROOT.gStyle.SetOptStat(0)

    c = ROOT.TCanvas(name,'',1200,1200)
    lMargin = 0.14
    rMargin = 0.10
    bMargin = 0.15
    tMargin = 0.10
    c.SetLeftMargin(lMargin)
    c.SetRightMargin(rMargin)
    c.SetTopMargin(tMargin)
    c.SetBottomMargin(bMargin)
    c.SetFrameBorderMode(0);
    c.SetBorderMode(0);
    c.SetBorderSize(0);
    return c

def doLegend(histos,labels,styles,corner="TR",textSize=0.035,legWidth=0.18,legBorder=False,nColumns=1):
    nentries = len(histos)
    (x1,y1,x2,y2) = (.85-legWidth, .7 - textSize*max(nentries-3,0), .90, .89)
    if corner == "TR":
        (x1,y1,x2,y2) = (.85-legWidth, .7 - textSize*max(nentries-3,0), .90, .89)
    elif corner == "TC":
        (x1,y1,x2,y2) = (.5, .7 - textSize*max(nentries-3,0), .5+legWidth, .89)
    elif corner == "TL":
        (x1,y1,x2,y2) = (.2, .7 - textSize*max(nentries-3,0), .2+legWidth, .89)
    elif corner == "BR":
        (x1,y1,x2,y2) = (.85-legWidth, .15 + textSize*max(nentries-3,0), .90, .25)
    elif corner == "BC":
        (x1,y1,x2,y2) = (.5, .15 + textSize*max(nentries-3,0), .5+legWidth, .25)
    elif corner == "BL":
        (x1,y1,x2,y2) = (.2, .23 + textSize*max(nentries-3,0), .33+legWidth, .35)
    leg = ROOT.TLegend(x1,y1,x2,y2)
    leg.SetNColumns(nColumns)
    leg.SetFillColor(0)
    leg.SetFillColorAlpha(0,0.6)  # should make the legend semitransparent (second number is 0 for fully transparent, 1 for full opaque)
    #leg.SetFillStyle(0) # transparent legend, so it will not cover plots (markers of legend entries will cover it unless one changes the histogram FillStyle, but this has other effects on color, so better not touching the FillStyle)
    leg.SetShadowColor(0)
    if not legBorder:
        leg.SetLineColor(0)
        leg.SetBorderSize(0)  # remove border  (otherwise it is drawn with a white line, visible if it overlaps with plots
    leg.SetTextFont(42)
    leg.SetTextSize(textSize)
    for (plot,label,style) in zip(histos,labels,styles): leg.AddEntry(plot,label,style)
    leg.Draw()
    ## assign it to a global variable so it's not deleted
    global legend_
    legend_ = leg
    return leg


class GBRLikelihoodTrainer:
    def set_defaults(self,params):
        self.tree_name = params['tree_name']
        self.target = params['target']
        self.var = params['inputs']
        self.cuts_base = params['selection']
        self.X_test = None
        self.y_test = None
        self.training = False
        
    def __init__(self,params):
        self.set_defaults(params)
        training_keys = ['n_estimators','max_depth','min_samples_split','learning_rate']
        self.training_params = {k: params[k] for k in training_keys}

    def name(self):
        return "{args.base_name}_{args.vars_name}_{args.cuts_name}_vgem{args.vgem1}V_ntrees{args.ntrees}".format(args=self)

    def variables(self):
        return self.var.split("|")
    
    def get_dataset(self,rfile):
        tfile = ROOT.TFile.Open(rfile)
        tree = tfile.Get(self.tree_name)

        # so target is always the first variable
        variables = [self.target] + self.var.split("|")
        print("List of variables = ",variables)
        dataset = tree2array(tree,variables,object_selection={self.cuts_base : variables})
        tfile.Close()
        
        conc = np.stack(dataset[0],axis=-1)
        for i in range(1,len(dataset)):
            conc = np.append(conc,np.stack(dataset[i],axis=-1),axis=0)
        X = conc[:,1:]
        y = conc[:,0]
        return X,y

    def train_model(self,X,y,options):
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=13)

        self.models_ = {}

        # MEAN SQUARE ERRORS REGRESSION
        print("===> Training mean square errors regression...")
        reg_ls = ensemble.GradientBoostingRegressor(loss='squared_error',
                                                    **self.training_params)
        self.models_["mse"] = reg_ls.fit(X_train, y_train)
        mse = mean_squared_error(y_test, reg_ls.predict(X_test))
        print("The mean squared error (MSE) on test set: {:.4f}".format(mse))

        # QUANTILE REGRESSION
        print("===> Now training quantiles regression...")
        alphas = [] if options.cvOnly==True else [0.05, 0.5, 0.95]
        for alpha in alphas:
            print("\t### Quantile = ",alpha)
            reg = ensemble.GradientBoostingRegressor(loss='quantile', alpha=alpha,
                                                     **self.training_params)
            self.models_["q%1.2f" % alpha] = reg.fit(X_train, y_train)
        
        self.X_test = X_test
        self.y_test = y_test
        self.training = True

    def plot_training(self):
        test_score = np.zeros((self.training_params['n_estimators'],), dtype=np.float64)
        for i, y_pred in enumerate(self.models_["mse"].staged_predict(self.X_test)):
            test_score[i] = self.models_["mse"].loss_(self.y_test, y_pred)

        fig = plt.figure(figsize=(6, 6))
        plt.subplot(1, 1, 1)
        plt.title('Deviance')
        plt.plot(np.arange(self.training_params['n_estimators']) + 1, self.models_["mse"].train_score_, 'b-',
                 label='Training Set Deviance')
        plt.plot(np.arange(self.training_params['n_estimators']) + 1, test_score, 'r-',
                 label='Test Set Deviance')
        plt.legend(loc='upper right')
        plt.xlabel('Boosting Iterations')
        plt.ylabel('Deviance')
        plt.savefig('training.png',bbox_inches='tight', pad_inches=0)
      
        # plot variables importance
        feature_importance = self.models_["mse"].feature_importances_
        sorted_idx = np.argsort(feature_importance)
        pos = np.arange(sorted_idx.shape[0]) + .5
        fig = plt.figure(figsize=(12, 6))
        plt.subplot(1, 2, 1)
        plt.barh(pos, feature_importance[sorted_idx], align='center')
        plt.yticks(pos, np.array(self.variables())[sorted_idx])
        plt.title('Feature Importance (MDI)')
        
        result = permutation_importance(self.models_["mse"], self.X_test, self.y_test, n_repeats=10,
                                        random_state=42, n_jobs=2)
        sorted_idx = result.importances_mean.argsort()
        plt.subplot(1, 2, 2)
        plt.boxplot(result.importances[sorted_idx].T,
                    vert=False, labels=np.array(self.variables())[sorted_idx])
        plt.title("Permutation Importance (test set)")
        plt.savefig('variables_importance.png')

        
    def save_models(self,prefix):
        for k in self.models_:
            filename = "{prefix}_{k}.sav".format(prefix=prefix.split('.')[0].replace(" ",""),k=k)
            print ("Saving model {k} in file {f}".format(k=k,f=filename))
            joblib.dump(self.models_[k],filename)
        print("Models saved.")

    def get_testdata(self):
        return self.X_test, self.y_test

    def test_models(self,recofile,prefix):
        print ("Test the saved models on the input file: ",recofile)
        # use the ones of the object if testing after training (to ensure orthogonality wrt training sample)
        if self.training:
            X_test = self.X_test
            y_test = self.y_test
        else:
            X,y = self.get_dataset(recofile)
            X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.95, random_state=13)

        hist = ROOT.TH1F('hist','',50,0.,2.0)
        hist.GetXaxis().SetTitle('E/E^{peak}_{raw}')
        hist.GetYaxis().SetTitle('Events')
        
        c = getCanvas('c')
        hists = {}
        maxy=-1
        for k in ['mse','q0.50']:
            filename = "{prefix}_{k}.sav".format(prefix=prefix.split('.')[0].replace(" ",""),k=k)
            print("Sanity check of the saved GBR model in the output file ",filename)
            model = joblib.load(filename)
            result = model.score(X_test, y_test)
            print("The precision on the test data is: ",result)
            y_pred = model.predict(X_test)
            hists[k] = hist.Clone('hist_{k}'.format(k=k))
            #print ("Regressed values = ",y_pred)
            fill_hist(hists[k],y_pred)
            maxy = max(maxy,hists[k].GetMaximum())
            
        hists["uncorr"] = hist.Clone('hist_uncorr')
        fill_hist(hists["uncorr"],y_test)
        labels = {'uncorr': "raw ({rms:1.2f}%)".format(rms=hists['uncorr'].GetRMS()),
                  'mse': 'regr. mean ({rms:1.2f}%)'.format(rms=hists['mse'].GetRMS()),
                  'q0.50': 'regr. median ({rms:1.2f}%)'.format(rms=hists['q0.50'].GetRMS())
                  }

        colors = {'uncorr': ROOT.kRed, 'mse': ROOT.kCyan, 'q0.50': ROOT.kBlack}
        styles = {'uncorr': 3005, 'mse': 3004, 'q0.50': 0}
        arr_hists = []; arr_styles = []; arr_labels = []
        for i,k in enumerate(colors):
            drawopt = '' if i==0 else 'same'
            hists[k].SetLineColor(colors[k])
            hists[k].SetFillColor(colors[k])
            hists[k].SetFillStyle(styles[k])
            hists[k].SetLineWidth(2)
            hists[k].SetMaximum(1.5 * maxy)
            hists[k].Draw("hist {opt}".format(opt=drawopt))
            # for the legend
            arr_hists.append(hists[k])
            arr_labels.append(labels[k])
            arr_styles.append('l')
        legend = doLegend(arr_hists,arr_labels,arr_styles,corner="TL")
        c.SaveAs("energy.png")
        
        
if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser(usage='%prog input.root params_gbrtrain.txt [opts] ')
    parser.add_option('-o', '--outname', dest='outname', default='gbrLikelihood', type='string', help='prefix for the output name of the regression models')
    parser.add_option('-a', '--apply-only', dest='applyOnly', action='store_true', default=False,  help='only apply regression on saved models and the data of the input file. Do not train')
    parser.add_option('--cv', '--central-value-only', dest='cvOnly', action='store_true', default=False,  help='Train only the central value regression, not the uncertainties.')
    (options, args) = parser.parse_args()

    recofile = args[0]
    config = open(args[1], "r")
    params = eval(config.read())

    GBR = GBRLikelihoodTrainer(params)
    if options.applyOnly == False:
        X,y = GBR.get_dataset(recofile)
        print("Dataset loaded from file ",args[0], " Now train the model.")
    
        GBR.train_model(X,y,options)
        print("GBR likelihood computed. Now plot results and control plots")
    
        GBR.plot_training()
        GBR.save_models(options.outname)

    # now test the model (this is to test that the model was saved correctly)
    GBR.test_models(recofile,options.outname)
    
    print("DONE.")
    
    

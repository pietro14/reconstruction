# USAGE training all regressions: python gbr_trainer.py ../prod-winter23-20230214/Merged/merged_feruns_8882_9857.root params_gbrtrain.txt -f ../prod-winter23-20230214/Merged/merged_feruns_8882_9857_Friend.root --savePanda regrdata_lngs_run2.pk
# USAGE testing only: python gbr_trainer.py ../prod-winter23-20230214/Merged/merged_feruns_8882_9857.root params_gbrtrain.txt -f ../prod-winter23-20230214/Merged/merged_feruns_8882_9857_Friend.root --loadPanda regrdata_lngs_run2.pkl -a

import ROOT
ROOT.gROOT.SetBatch(True)

import numpy as np
import awkward as ak
import uproot
import pandas as pd

import matplotlib.pyplot as plt
import numpy as np
import joblib

from sklearn import ensemble
from sklearn.inspection import permutation_importance
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split

def fill_hist(hist,arr):
    for i in range(len(arr)):
        hist.Fill(arr[i])

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
        self.target = params['target']
        self.tree_name = params['tree_name']
        self.var = params['inputs']
        self.varfriend = params['inputsfriend']
        self.regr_vars = params['regr_vars']
        self.X_test = None
        self.y_test = None
        self.training = False
        self.verbose = params['verbose']
        
    def __init__(self,paramsfile):
        config = open(paramsfile, "r")
        params = eval(config.read())
        self.set_defaults(params)
        training_keys = ['n_estimators','max_depth','min_samples_split','learning_rate']
        self.training_params = {k: params[k] for k in training_keys}

    def name(self):
        return "{args.base_name}_{args.vars_name}_{args.cuts_name}_vgem{args.vgem1}V_ntrees{args.ntrees}".format(args=self)

    def variables(self):
        return self.regr_vars.split("|")
    
    def get_dataset(self,rfile,friendrfile=None,firstEvent=None,lastEvent=None,savePanda=None,loadPanda=None,addCuts={}):
        variables_events = self.var.split("|")
        variables_friends = self.varfriend.split("|")
        # add regression inputs, only the variables which are function of the others and remove the eventual duplicates
        regr_inputs = self.regr_vars.split("|")
        variables_events = list(set(variables_events + regr_inputs))

        if not loadPanda:
            print ("Loading events from file %s and converting to numpy arrays for training. It might take time..." % rfile)
            events = uproot.open(rfile)
            if friendrfile:
                friends = uproot.open(friendrfile)
            if self.verbose: print ("---> Now loading main tree %s..." % rfile)
            data_main = events[self.tree_name].arrays(variables_events,library="pd",entry_start=firstEvent,entry_stop=lastEvent)
            if self.verbose: print ("---> Now loading friend tree %s ..." % friendrfile)
            data_friend = friends["Friends"].arrays(variables_friends,library="pd",entry_start=firstEvent,entry_stop=lastEvent)
            if len(data_main.index)!=len(data_friend.index): RuntimeError("Number of entries in the main tree = %d and in the friend tree = %d don't match " %(len(data_main.index),len(data_friend.index)))
            if self.verbose: print ("---> Now attaching main and friend pandas...")
            data = pd.concat([data_main,data_friend],axis=1) # Panda dataframe
            if self.verbose > 0:
                print (" ~~~~~~~~~~~~~~ DATA BEFORE SELECTION ~~~~~~~~~~~~~~~~~~ ")
                print (data)
                print (" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")
     
            if self.verbose: print ("---> Now calculating target variable and attaching to panda...")
            if self.target=='sc_trueint':
                data["target"] = data.apply(lambda row: row.sc_trueint, axis=1)
            elif self.target=='sc_truez':
                data["target"] = data.apply(lambda row: row.sc_truez, axis=1)
            else:
                RuntimeError("Target variable %s not foreseen. Either pass 'sc_trueint' or 'sc_truez'. Exiting." % self.target)
                
            if savePanda:
                if self.verbose: print ("---> Now saving selected clusters to panda into pikle file %s..." % savePanda)
                data.to_pickle(savePanda)
        else:
            data = pd.read_pickle(loadPanda)

        ### hardcoded, move to configuration
        if self.verbose: print ("---> Now applying selection to panda...")
        data_sel = data[(data['sc_trueint']>0)&(data['sc_integral']>1500)&(data['sc_rms']>8)&(data['sc_tgausssigma']*0.152>0.3)&(np.hypot(data['sc_xmean']-2304/2,data['sc_ymean']-2304/2)<900)]
        #if self.target=='sc_truez':
        #    data_sel = data_sel[(data_sel['sc_hv']==440)]
        if len(addCuts):
            for k,v in addCuts.items():
                print ("Adding selection: %d < %s <= %d " % (v[0],k,v[1]))
                data_sel = data_sel[(data_sel[k]>=v[0])&(data_sel[k]<v[1])]
        if self.verbose > 0:
            print (" ~~~~~~~~~~~~~~ DATA AFTER SELECTION ~~~~~~~~~~~~~~~~~ ")
            print (data_sel)
            print (" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")
            
        if self.verbose:
            print("List of regression variables = ",regr_inputs)
            print ("---> Now select only regression variables...")
        data_regr = data_sel[regr_inputs+["target"]]     

        if self.verbose > 0:
            print (" ~~~~~~~~~~~~~~ DATA REGRESSION ~~~~~~~~~~~~~~~~~ ")
            print (data_regr)
            print (" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")

        X = data_regr.to_numpy()[:,:-1]
        y = data_regr.to_numpy()[:,-1]
        self.rawyindex = data_regr.columns.get_loc("sc_integral") if self.target=="sc_trueint" else -1 # for Z there is not a raw estimate
        return X,y

    def train_model(self,X,y,options):
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.5, random_state=13)

        self.models_ = {}

        # MEAN SQUARE ERRORS REGRESSION
        print("===> Training mean square errors regression...")
        reg_ls = ensemble.GradientBoostingRegressor(loss='ls',
                                                    **self.training_params)
        self.models_["mse"] = reg_ls.fit(X_train, y_train)
        mse = mean_squared_error(y_test, reg_ls.predict(X_test))
        print("The mean squared error (MSE) on test set: {:.4f}".format(mse))

        # QUANTILE REGRESSION
        if not options.cvOnly:
            print("===> Now training quantiles regression...")
            alphas = [0.05, 0.5, 0.95]
        
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
        for ext in ['pdf','png']:
            plt.savefig('training_%s.%s'%(self.target.replace('sc_',''),ext),bbox_inches='tight', pad_inches=0)
      
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
        for ext in ['pdf','png']:
            plt.savefig('variables_importance_%s.%s'%(self.target.replace('sc_',''),ext))

        
    def save_models(self,prefix):
        for k in self.models_:
            filename = "{prefix}_{target}_{k}.sav".format(prefix=prefix.split('.')[0].replace(" ",""),target=self.target.replace('sc_',''),k=k)
            print ("Saving model {k} in file {f}".format(k=k,f=filename))
            joblib.dump(self.models_[k],filename)
        print("Models saved.")

    def get_testdata(self):
        return self.X_test, self.y_test

    def test_models(self,recofile,options,panda=None):
        prefix=options.outname
        print ("Test the saved models on the input file: ",recofile)
        # use the ones of the object if testing after training (to ensure orthogonality wrt training sample)
        if self.training:
            X_test = self.X_test
            y_test = self.y_test
        else:
            X,y = self.get_dataset(recofile,options.friend,loadPanda=panda)
            X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.95, random_state=13)

        xmin,xmax=(0,2) if self.target=="sc_trueint" else (-25,25)
        hist = ROOT.TH1F('hist','',50,xmin,xmax)        
        if self.target=="sc_trueint":
            hist.GetXaxis().SetTitle('E / E_{true}')
        else:
            hist.GetXaxis().SetTitle('Z - Z_{true}')
        hist.GetYaxis().SetTitle('Events')
        
        c = getCanvas('c')
        hists = {}
        maxy=-1
        alphas = ['mse'] if options.cvOnly==True else ['mse', 'q0.05', 'q0.50', 'q0.95']
        for k in alphas:
            filename = "{prefix}_{target}_{k}.sav".format(prefix=prefix.split('.')[0].replace(" ",""),target=self.target.replace('sc_',''),k=k)
            print("Sanity check of the saved GBR model in the output file ",filename)
            model = joblib.load(filename)
            result = model.score(X_test, y_test)
            print("The precision on the test data is: ",result)
            y_pred = model.predict(X_test)
            hists[k] = hist.Clone('hist_{k}'.format(k=k))
            print ("True values = ",y_test)
            print ("Raw values = ",X_test[:,self.rawyindex])
            print ("Regressed values = ",y_pred)
            if self.target=='sc_trueint':
                YoYt = np.divide(y_pred,y_test) # Predicted Y over Ytrue
                print ("Pred/True = ",YoYt)
            else:
                YoYt = np.subtract(y_pred,y_test) # Predicted Y - Ytrue
                print ("Pred-True = ",YoYt)                
            fill_hist(hists[k],YoYt)
            maxy = max(maxy,hists[k].GetMaximum())

        # for z there is not a raw estimate
        if self.target=='sc_trueint':
            hists["uncorr"] = hist.Clone('hist_uncorr')
            YrawoYt = np.divide(X_test[:,self.rawyindex],y_test)
            print ("Raw/True = ",YrawoYt)
            fill_hist(hists["uncorr"],YrawoYt)
            
        labels = {'mse': 'regr. mean ({rms:1.2f}%)'.format(rms=hists['mse'].GetRMS())}
        if 'uncorr' in hists: labels['uncorr'] = "raw ({rms:1.2f}%)".format(rms=hists['uncorr'].GetRMS())
        if 'q0.50' in hists: labels['q0.50'] = 'regr. median ({rms:1.2f}%)'.format(rms=hists['q0.50'].GetRMS())
        colors = {'uncorr': ROOT.kRed, 'mse': ROOT.kCyan, 'q0.50': ROOT.kBlack}
        styles = {'uncorr': 3005, 'mse': 3004, 'q0.50': 0}
        arr_hists = []; arr_styles = []; arr_labels = []
        for i,k in enumerate(hists):
            drawopt = '' if i==0 else 'same'
            if k not in colors.keys(): continue 
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
        for ext in ['pdf','png']: 
            c.SaveAs("%s.%s"%(self.target.replace('sc_',''),ext))
        
        
if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser(usage='%prog input.root params_gbrtrain.txt [opts] [-f friend.root] ')
    parser.add_option('-o',   '--outname', dest='outname', default='gbrLikelihood', type='string', help='prefix for the output name of the regression models')
    parser.add_option('-a',   '--apply-only', dest='applyOnly', action='store_true', default=False,  help='only apply regression on saved models and the data of the input file. Do not train')
    parser.add_option('--cv', '--central-value-only', dest='cvOnly', action='store_true', default=False,  help='Train only the central value regression, not the uncertainties.')
    parser.add_option('-f',   '--friend', dest='friend', default=None, type='string', help='file to be used as friend (eg for Etrue)')
    parser.add_option(        '--savePanda', dest='savePanda', default=None, type='string', help='file where to store the regression data as panda dataframe for re-use')
    parser.add_option(        '--loadPanda', dest='loadPanda', default=None, type='string', help='file with regression data as panda dataframe to be loaded instead of the full ROOT files')
    (options, args) = parser.parse_args()

    recofile = args[0]
    paramsfile = args[1]

    GBR = GBRLikelihoodTrainer(paramsfile)
    if options.applyOnly == False:
        X,y = GBR.get_dataset(recofile,friendrfile=options.friend,savePanda=options.savePanda,loadPanda=options.loadPanda)
        print("Dataset loaded from file ",args[0], " Now train the model.")

        GBR.train_model(X,y,options)
        print("GBR likelihood computed. Now plot results and control plots")
        
        GBR.plot_training()
        GBR.save_models(options.outname)

    # now test the model (this is to test that the model was saved correctly)
    GBR.test_models(recofile,options,panda=options.loadPanda)
    
    print("DONE.")
    
    

from .tfTool import TFTool

class finalDNN:
    def __init__(self):
        self.outVars = []
        
        varorder = ["density","length","slimness","lgsigma","tgsigma","latrms","longrms","size","nhits"]
        cats = ['pred_nr','pred_er','pred_other']
        self._MVA = TFTool('DNN','data/mvas/trained_model_A_9vars_DATA_31-03-2021.h5',
                           self.getVars(), cats, varorder)
        self.outVars.extend( ['DNN%s_' % x for x in cats])

    def getVars(self):
        return {     "density"          : lambda ev,isc : ev.sc_integral[isc]/ev.sc_nhits[isc] if  ev.sc_integral[isc]>0 else -9, # 0
                     "length"           : lambda ev,isc : ev.sc_length[isc],                                                      # 1
                     "slimness"         : lambda ev,isc : ev.sc_width[isc]/ev.sc_length[isc],                                     # 2
                     "lgsigma"          : lambda ev,isc : ev.sc_lgausssigma[isc],                                                 # 3
                     "tgsigma"          : lambda ev,isc : ev.sc_tgausssigma[isc],                                                 # 4
                     "latrms"           : lambda ev,isc : ev.sc_latrms[isc],                                                      # 5
                     "longrms"          : lambda ev,isc : ev.sc_longrms[isc],                                                     # 6
                     "size"             : lambda ev,isc : ev.sc_size[isc],                                                        # 7 
                     "nhits"            : lambda ev,isc : ev.sc_nhits[isc],                                                       # 8
                }

    def analyze(self,event,isc):
        ret = {}
        worker = self._MVA
        for node,val in worker(event,isc).items():
            ret[node] = val
        return ret



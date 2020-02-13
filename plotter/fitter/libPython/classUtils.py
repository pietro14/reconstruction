def mkdir(directory):
    if not os.path.isdir(directory):
        os.makedirs(directory) 

class fitSample:
    def __init__(self, sName, path, cut = None, nEvts = -1 ):
        self.path = path
        self.name = sName
        self.cut     = cut
        self.nEvts   = nEvts
        self.fitOutput = 'fitresults.root'
        
    def clone(self):
        return copy.deepcopy(self)


    

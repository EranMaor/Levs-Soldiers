import numpy as np
class Channel:
    def probsim(self, p):
        if np.random.random()<p:
            return True
        else:
            return False
    def __init__(self, p_ins, p_del, Dmin, Dmax, p_sub=np.eye(4)):    
        self.p_ins = p_ins
        self.p_del = p_del
        self.Dmin = Dmin
        self.Dmax = Dmax
        self.p_sub= p_sub
    def drift(self, x):
        d=0
        drifted_x=[]
        for ind in x:
            if d > self.Dmin:
                if  self.probsim(self.p_del):
                    d+=-1
                    continue
                elif d < self.Dmax and self.probsim((self.p_ins)/(1-self.p_del)):
                    d+=1
                    drifted_x.append(ind)
            elif self.probsim(self.p_ins):
                d+=1
                drifted_x.append(ind)
            drifted_x.append(ind)
        return(drifted_x)
    def encode(self, X):
        DX = [self.drift(x) for x in X]
        Y=[]
        for strand in DX:
            current=[]
            for ind in strand:
                nucleo = np.random.choice([0,1,2,3],p=self.p_sub[:, ind])
                current.append(nucleo)
            Y.append(current)
        return Y

        

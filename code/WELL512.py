"""
Alex Zolan
March 30, 2017
WELL512a random number generator

This is an implementation of the WELL512a RNG created by L'Ecuyer and 
Matsumoto, originally programmed in C.
"""

#!java -classpath ./ssj/lib/ssj.jar

import scipy
#import py4j
#from py4j.java_gateway import JavaGateway 


class WELL512(object):
    def __init__(self,state_file,W=32,R=16,P=0,M1=13,M2=9,M3=5,scen=0):
        self.state = scipy.zeros(R,dtype=scipy.uint32)
        self.W = W
        self.R = R
        self.P = P
        self.M1 = M1
        self.M2 = M2
        self.M3 = M3
        self.FACT = 2.32830643653869628906e-10
        states, state_is = self.readStates(state_file)
        self.state = states[scen]
        self.state_i = state_is[scen]
        
    def readStates(self,f):
        states = []
        state_is = []
        import csv
        reader = csv.reader(open(f,'rU'))
        for line in reader:
            state = []
            for i in range(16):
                state.append(scipy.uint32(line[i]))
            if len(line) > 16:
                state_i = int(line[-1])
            else: 
                state_i = 0
            states.append(state)
            state_is.append(state_i)
        return states, state_is
        
    def setState(self,idx):
        self.state = self.states[idx]
            
    def V0(self):
        return self.state_i
        
    def VM1(self):
        return (self.state_i + self.M1) & 0x0000000f 
           
    def VM2(self):
        return (self.state_i + self.M2) & 0x0000000f    
        
    def VM3(self):
        return (self.state_i + self.M3) & 0x0000000f
            
    def VRm1(self):
        return (self.state_i + 15) & 0x0000000f
    
    def VRm2(self):
        return (self.state_i + 14) & 0x0000000f
    
    def newV0(self):
        return (self.state_i + 15) & 0x0000000f
    
    def newV1(self):
        return self.state_i
        
    def newVRm1(self):
        return (self.state_i + 14) & 0x0000000f
    
    def MAT0POS(self,t,v): 
        return scipy.uint32(v^(v>>t))
        
    def MAT0NEG(self,t,v): 
        return scipy.uint32(v^(v<<(-(t))))
        
    def MAT3NEG(self,t,v): 
        return scipy.uint32(v<<(-(t)))
        
    def MAT4NEG(self,t,b,v): 
        return scipy.uint32(v ^ ((v<<(-(t))) & b))
        
    def getVariate(self):
        z0 = scipy.uint32(self.state[self.VRm1()])
        z1 = scipy.uint32(self.MAT0NEG(-16, self.V0())) ^ scipy.uint32(self.MAT0NEG(-15, self.state[self.VM1()]))
        z2 = scipy.uint32(self.MAT0POS(11, self.state[self.VM2()]))
        self.state[self.newV1()] = z1 ^ z2
        self.state[self.newV0()] = self.MAT0NEG(-2, z0) ^ self.MAT0NEG(-18, z1) ^ self.MAT3NEG(-28, z2) ^ self.MAT4NEG(-5, scipy.uint32(0xda442d24), self.state[self.newV1()])
        self.state_i = (self.state_i + 15) & 0x0000000f
        return self.state[self.state_i] * self.FACT
        
    def getVariates(self,num_samples):
        arr = scipy.zeros(num_samples,dtype=scipy.float64)
        for i in range(num_samples):
            arr[i] = self.getVariate()
        return arr
        
    def getPermutation(self,num_items):
        arr = self.getVariates(num_items)
        return arr.argsort()
        
    def GetState(self):
        return self.state, self.state_i
        
    
    
if __name__ == "__main__":
    w = WELL512("rngstates.csv")
#    print w.state
#    w.setState(88)
#    print w.state   
    
        
        
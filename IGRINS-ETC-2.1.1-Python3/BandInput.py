'''
Created on Aug  24, 2011

@Author: Jongmin LEE and Huynh Anh Nguyen LE
'''
from tkinter import *

class CBand():
    '''
    classdocs
    '''
    

    def __init__(self, master, rownum, colnum):
        '''
        Constructor
        '''
        self.master = master
        self.rownum = rownum
        self.colnum = colnum
             
        self.BANDS = [self.e_h, self.e_k] = [Entry(self.master, width=7, justify=RIGHT),
                                   Entry(self.master, width=7, justify=RIGHT)]
        
        self.VALUE = [e1, e2] = [StringVar(), StringVar()]
              
    def ShowEntry(self, mode):
        for i, m in enumerate(self.BANDS):
            if mode == "active" :                
                m.config(state="normal")
        
            elif mode == "disable" : 
                m.config(state="disable")
            
            elif mode == "readonly" : 
                m.config(state="readonly")
            
            m.grid(row=self.rownum, column=self.colnum+i, padx=1, pady=3)
    
    def GetEntry(self):
        for i, m in enumerate(self.BANDS):
            self.VALUE[i] = m.get()
        return float(self.VALUE[0]), float(self.VALUE[1])
    
    def ShowValue(self, v):
        for i, m in enumerate(self.BANDS):
            m.delete(0, END)
            m.insert(0, v[i])
            
    def GetActiveEntryValue(self, type):
       
        if type == "h" : 
            value = self.e_h.get()
            
        elif type == "k" : 
            value = self.e_k.get()
            
        return value
        
            
    def ActiveEntry(self, type):
        for m in self.BANDS:
            m.config(state="disable")
        
        #if type == "h" : 
        self.e_h.config(state="normal")
            
        #elif type == "k" : 
        self.e_k.config(state="normal")
            
        
class CBandLabel(CBand):
    '''
    classdocs
    '''
    

    def __init__(self, master, rownum, colnum):
        '''
        Constructor
        '''
        MCBand.__init__(self, master, rownum, colnum)
        
        self.BANDSLABEL = [self.t_h, self.t_k] = [Label(self.master, text="H = ", width=4, anchor="center"),
                                        Label(self.master, text="K = ", width=4, anchor="center")]
              
    def ShowEntry(self, mode):
        MCBand.ShowEntry(self, mode)
        
        for i, n in enumerate(self.BANDSLABEL):
            if mode == "active" :                
                n.config(state="active")
        
            elif mode == "disable" : 
                n.config(state="disable")
            
            else :
                pass
                
            n.grid(row=self.rownum+i, column=self.colnum-1, padx=1)
    
    #def GetEntry(self):
    #    v1, v2 = MCBand.GetEntry(self)
    #    return v1, v2
    
    def ShowValue(self, v):
        MCBand.ShowValue(self, v)
        return

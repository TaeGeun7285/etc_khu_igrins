"""
Created on Sep 19, 2011

@Author: Jongmin LEE and Huynh Anh Nguyen LE
"""

from tkinter import *
from Display import *
from Constants import *
    
class CFrame(Frame):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''   
        Frame.__init__(self)
#        self.master.title("IGRINS Exposure Time Calculator")
#        self.master.iconname("IGRINS Exposure Time Calculator")
        self.master.title(ETC_version)
        self.master.iconname(ETC_version)

        
        self.Main()    
        self.Draw()
                   
    def Main(self):
        mainframe = PanedWindow(self.master, orient="vertical")
        self.mode = LabelFrame(mainframe, text="Calculation Mode Selection")
        self.parameter = LabelFrame(mainframe, text="User Input Parameters")
        self.resultsframe = Frame(mainframe)
                
        mainframe.add(self.mode)
        mainframe.add(self.parameter)
        mainframe.add(self.resultsframe)
        mainframe.grid(padx=30, pady=30)
        
        self.Results()
                
    def Results(self):   
        mainframe = PanedWindow(self.resultsframe, orient="horizontal")
        self.button = Frame(mainframe)
        self.results = LabelFrame(mainframe, text="Results")
                
        mainframe.add(self.button)
        mainframe.add(self.results)
        mainframe.grid()
        
    def Draw(self):   
        self.draw = CDisplay(self.mode, self.parameter, self.button, self.results)

        self.draw.Mode()
        self.draw.Parameter()
        self.draw.Button()
        self.draw.Results()

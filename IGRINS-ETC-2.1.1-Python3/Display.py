"""
Created on Sep 19, 2011

@Author: Jongmin LEE and Huynh Anh Nguyen LE
"""

from BandInput import *
from Functions import *
from tkinter import *

from Init_Value import *

class CDisplay():
    '''
    classdocs
    '''
    

    def __init__(self, mode, parameter, button, results):
        '''
        Constructor
        '''
        self.mode = StringVar()
        self.mode.set("single")
        self.mag = StringVar()
        
        self.modframe = mode
        self.parframe = parameter
        self.butframe = button
        #self.resframe = results
        self.x, self.y = (StringVar(), StringVar())
        
        #self.SINGLE = [(self.t_point, self.magband), 
        #               (self.t_snr, self.snrband)] = [(Label(self.parframe, width=4, text = "Mag."), CBandLabel(self.parframe, 5, 5)),
        #                                              (Label(self.resframe, width=4, text = "S/N"), CBandLabel(self.resframe, 2, 4))]
        

        self.SINGLE = (self.t_mag_h, self.e_mag_h, self.t_mag_k,
                     self.e_mag_k) = (Label(self.parframe, width=20, text = "H [mag]     = ", anchor=CENTER, justify=CENTER), \
                                      Label(self.parframe, width=7, anchor=CENTER, justify=CENTER, textvariable=self.x),
                                        Label(self.parframe, width=20, text = "K [mag]     = ", anchor=CENTER, justify=CENTER), \
                                      Label(self.parframe, width=7, anchor=CENTER, justify=CENTER, textvariable=self.y))


        self.PLOT = (self.t_mag_min, self.e_mag_min, self.t_mag_max,
                     self.e_mag_max) = (Label(self.parframe, width=1, text = ""), Entry(self.parframe, width=7, justify=RIGHT, textvariable=DoubleVar(value=minimum_)),
                                        Label(self.parframe, width=1, text = ""), Entry(self.parframe, width=7, justify=RIGHT, textvariable=DoubleVar(value=maximum_)))
    
        #self.MAGNITUDE = [(self.t_mag, self.signoiseband)] = [(Label(self.resframe, width=8, text = "Mag."), CBandLabel(self.resframe, 2, 7))]    
        
        self.wavemin = CBand(self.parframe, 6, 9) 
        self.wavemax = CBand(self.parframe, 8, 9)
                                           
    def Mode(self):
        Label(self.modframe, width=35, text = "").grid(row=1, column=1, columnspan=1, sticky=W)
        
        single = Radiobutton(self.modframe, text="S/N Calculation", variable=self.mode, value="single", command=self.ModeSelect)
        single.grid(row=1, column=2, padx=1, sticky=W)
        
        plot = Radiobutton(self.modframe, text="S/N vs. Magnitude", variable=self.mode, value="plot", command=self.ModeSelect)
        plot.grid(row=1, column=4, padx=1, sticky=W)
        
        plot = Radiobutton(self.modframe, text="S/N vs. Wavelength", variable=self.mode, value="plot_wave", command=self.ModeSelect)
        plot.grid(row=2, column=4, padx=1, sticky=W)
        
        plot = Radiobutton(self.modframe, text="Signal with Noise vs. Wavelength", variable=self.mode, value="plot_signal", command=self.ModeSelect)
        plot.grid(row=3, column=4, padx=1, sticky=W)

        plot = Radiobutton(self.modframe, text="Noise vs. Wavelength", variable=self.mode, value="plot_noise", command=self.ModeSelect)
        plot.grid(row=4, column=4, padx=1, sticky=W)

        Label(self.modframe, width=18, text = "").grid(row=3, column=3, columnspan=1, sticky=W)
        
        #magnitude = Radiobutton(self.modframe, text="Limiting Magnitude", variable=self.mode, value="magnitude", command=self.ModeSelect)
        #magnitude.grid(row=1, column=4, padx=1, sticky=W)
                	
    def Parameter(self):
        self.t_pwv = Label(self.parframe, width=20, text="PWV [mm]")
        self.t_pwv.grid(row=1, column=1, columnspan=1, sticky=W)
        self.e_pwv = Entry(self.parframe, width=7, justify=RIGHT, textvariable=DoubleVar(value=pwv_))
        self.e_pwv.grid(row=1, column=2, padx=1, pady=1, columnspan=1, sticky=W)

		
        Label(self.parframe, width=20, text="Exp.T_ [sec]").grid(row=2, column=1, columnspan=1, sticky=W)
        self.exptime = Entry(self.parframe, width=7, justify=RIGHT, textvariable=DoubleVar(value=exptime_))
        self.exptime.grid(row=2, column=2, padx=1, pady=1, columnspan=1, sticky=W)
        
        Label(self.parframe, width=20, text="Num of Exp.").grid(row=3, column=1, columnspan=1, sticky=W)
        self.expnumber = Entry(self.parframe, width=7, justify=RIGHT, textvariable=DoubleVar(value=expnumber_))
        self.expnumber.grid(row=3, column=2, padx=1, pady=1, columnspan=1,  sticky=W)

        self.t_kmag = Label(self.parframe, width=20, text="K [mag]", anchor=CENTER, justify=CENTER)
        self.t_kmag.grid(row=1, column=4, columnspan=1, sticky=W)
        self.e_kmag = Entry(self.parframe, width=7, justify=RIGHT, textvariable=DoubleVar(value=kmag_))
        self.e_kmag.grid(row=1, column=5, padx=1, pady=1, columnspan=1, sticky=W)

        
        self.t_temperature = Label(self.parframe, width=20, text="T_bb [K]", anchor=CENTER, justify=CENTER)
        self.t_temperature.grid(row=2, column=4, columnspan=1, sticky=W)
        self.e_temperature = Entry(self.parframe, width=7, justify=RIGHT, textvariable=DoubleVar(value=temperature_))
        self.e_temperature.grid(row=2, column=5, padx=1, pady=1, columnspan=1, sticky=W)


        
        Label(self.parframe, width=5, text = " ").grid(row=1, column=9, columnspan=1, sticky=W)
        
        
        #self.plot_H = Radiobutton(self.parframe, text="", variable=self.mag, value="plot_wave_H", command=self.MagSelect)
        #self.plot_H.grid(row=5, column=8, padx=1, sticky=W)
        #self.plot_K = Radiobutton(self.parframe, text="", variable=self.mag, value="plot_wave_K", command=self.MagSelect)
        #self.plot_K.grid(row=6, column=8, padx=1, sticky=W)


        self.plot_H = Radiobutton(self.parframe, text="", variable=self.mag, value="plot_wave_H", command=self.MagSelect)
        self.plot_H.grid(row=6, column=8, padx=1, sticky=W)
        self.plot_K = Radiobutton(self.parframe, text="", variable=self.mag, value="plot_wave_K", command=self.MagSelect)
        self.plot_K.grid(row=8, column=8, padx=1, sticky=W)
        Label(self.parframe, width=7, text = "[1.4   - ").grid(row=5, column=9, columnspan=1, sticky=W)
        Label(self.parframe, width=5, text = "1.9]").grid(row=5, column=10, columnspan=1, sticky=W)
        Label(self.parframe, width=7, text = "[1.9   -  ").grid(row=7, column=9, columnspan=1, sticky=W)
        Label(self.parframe, width=5, text = "2.5]").grid(row=7, column=10, columnspan=1, sticky=W)
               
        Label(self.parframe, width=1, text = " ").grid(row=1, column=3, columnspan=1, sticky=W)

        #self.t_point.grid(row=4, column=5)
        #self.magband.ShowEntry("enable")  

        self.t_mag_h.grid(row=5, column=4, columnspan=1, sticky=W)        
        self.e_mag_h.grid(row=5, column=5, sticky=W)

        self.t_mag_k.grid(row=6, column=4, columnspan=1, sticky=W)                
        self.e_mag_k.grid(row=6, column=5, sticky=W)

               
        self.t_mag_min.grid(row=2, column=9, columnspan=1, sticky=W)        
        self.e_mag_min.grid(row=3, column=9, sticky=W)

        self.t_mag_max.grid(row=2, column=10, columnspan=1, sticky=W)                
        self.e_mag_max.grid(row=3, column=10, sticky=W)


        self.t_waverange = Label(self.parframe, width=20, text="Wavelength Range [um]")
        self.t_waverange.grid(row=4, column=8, columnspan=7, sticky=W)
        #self.t_waverange = Label(self.parframe, width=15, text="(1.4 - 1.9) um", anchor=W, justify=LEFT)
    	#self.t_waverange.grid(row=5, column=11, columnspan=1, sticky=W)
        #self.t_waverange = Label(self.parframe, width=15, text="(1.9 - 2.5) um", anchor=W, justify=LEFT)
    	#self.t_waverange.grid(row=6, column=11, columnspan=1, sticky=W)

        self.t_magrange = Label(self.parframe, width=20, text="Mag. Range [mag]")
        self.t_magrange.grid(row=2, column=8, columnspan=5, sticky=W)
         
        
        #Label(self.parframe, width=3, text = "mag").grid(row=5, column=6, columnspan=1, sticky=W)
        #Label(self.parframe, width=3, text = "mag").grid(row=6, column=6, columnspan=1, sticky=W)

        Label(self.parframe, width=10, text = " ").grid(row=5, column=7, columnspan=1, sticky=W)
        Label(self.parframe, width=10, text = " ").grid(row=6, column=7, columnspan=1, sticky=W)

        self.wavemin.ShowEntry("enable")
        vlist = [waveminH_, waveminK_]
        self.wavemin.ShowValue(vlist)
        
        self.wavemax.ShowEntry("enable")
        vlist = [wavemaxH_, wavemaxK_]
        self.wavemax.ShowValue(vlist)  
                     
        self.t_emission = Label(self.parframe, width=20, text="Emission Line", anchor=W, justify=LEFT)
        self.t_emission.grid(row=7, column=5, columnspan=4, sticky=W)

        self.t_restwave = Label(self.parframe, width=20, text="Rest Wavelength [um]", anchor=W, justify=LEFT)
        self.t_restwave.grid(row=8, column=4, columnspan=1, sticky=W)
        self.e_restwave = Entry(self.parframe, width=7, justify=RIGHT, textvariable=DoubleVar(value=restwave_))
        self.e_restwave.grid(row=8, column=5, padx=1, pady=1, columnspan=1, sticky=W)
        
        self.t_lineflux = Label(self.parframe, width=20, text="Line Flux [1E-18 W m-2]", anchor=W, justify=LEFT)
        self.t_lineflux.grid(row=9, column=4, columnspan=1, sticky=W)
        self.e_lineflux = Entry(self.parframe, width=7, justify=RIGHT, textvariable=DoubleVar(value=lineflux_))
        self.e_lineflux.grid(row=9, column=5, padx=1, pady=1, columnspan=1, sticky=W)
        

        self.t_doppler = Label(self.parframe, width=20, text="Doppler Shift [Km s-1]", anchor=W, justify=LEFT)
        self.t_doppler.grid(row=10, column=4, columnspan=1, sticky=W)
        self.e_doppler = Entry(self.parframe, width=7, justify=RIGHT, textvariable=DoubleVar(value=doppler_))
        self.e_doppler.grid(row=10, column=5, padx=1, pady=1, columnspan=1, sticky=W)
        

        self.t_linewidth = Label(self.parframe, width=20, text="Line Width [Km s-1]", anchor=W, justify=LEFT)
        self.t_linewidth.grid(row=11, column=4, columnspan=1, sticky=W)
        self.e_linewidth = Entry(self.parframe, width=7, justify=RIGHT, textvariable=DoubleVar(value=linewidth_))
        self.e_linewidth.grid(row=11, column=5, padx=1, pady=1, columnspan=1, sticky=W)
        


    def Results(self):       
    #    Label(self.resframe, width=12, text = "").grid(row=1, column=2, columnspan=1, sticky=W)
    #    Label(self.resframe, width=18, text = "").grid(row=1, column=5, columnspan=1, sticky=W)
    #    Label(self.resframe, width=7, text = "").grid(row=1, column=8, columnspan=1, sticky=W)
        
        #self.t_snr.grid(row=1, column=4, rowspan=1, padx=1)
        #self.snrband.ShowEntry("readonly") 

        #self.t_mag.grid(row=1, column=7, rowspan=1, padx=1)
        #self.signoiseband.ShowEntry("readonly")
        
       self.ModeSelect() 
		
    def Button(self):
        Button(self.butframe, text="Run", width=8, command=self.Run).grid(row=1, column=1, padx=10, pady=15)
        Button(self.parframe, text="Template", width=8, command=self.Get_Mag).grid(row=3, column=5, padx=1, pady=3)
        
    def ModeSelect(self):
        if self.mode.get() == "single":
            #for t, b in self.MAGNITUDE :
            #    t.config(state='disable')
            #    b.ShowEntry("disable")
            
            for t in self.SINGLE :
                t.config(state='normal')
                

            #for t, b in self.SINGLE :
            #    t.config(state='active')
            #    b.ShowEntry("active")
            
            for f in self.PLOT :
                f.config(state='disable')
            
            self.t_pwv.config(state='disable')
            self.e_pwv.config(state='disable')
            #self.t_signalnoise.config(state='disable')
            #self.e_signalnoise.config(state='disable')
            self.t_waverange.config(state='disable')
            self.t_magrange.config(state='disable')
            self.t_kmag.config(state='normal')
            self.e_kmag.config(state='normal')
            self.t_temperature.config(state='normal')
            self.e_temperature.config(state='normal')
            #self.t_point.config(state='disable')
            
            self.wavemin.ShowEntry("disable")
            self.wavemax.ShowEntry("disable")    
			
            self.mag.set(" ")
            self.plot_H.config(state='disable')
            self.plot_K.config(state='disable')

            self.t_emission.config(state='disable')
            self.t_restwave.config(state='disable')
            self.e_restwave.config(state='disable')
            #self.u_restwave.config(state='disable')
            self.t_doppler.config(state='disable')
            self.e_doppler.config(state='disable')
            #self.u_doppler.config(state='disable')
            self.t_lineflux.config(state='disable')
            self.e_lineflux.config(state='disable')
            #self.u_lineflux.config(state='disable')
            self.t_linewidth.config(state='disable')
            self.e_linewidth.config(state='disable')
            #self.u_linewidth.config(state='disable')

        elif self.mode.get() == "plot":
            #for t, b in self.MAGNITUDE :
            #    t.config(state='disable')
            #    b.ShowEntry("disable")
			
            for t in self.SINGLE :
                t.config(state='disable')

            #for t, b in self.SINGLE :
            #    t.config(state='disable')
            #    b.ShowEntry("disable")
                
            for f in self.PLOT :
                f.config(state='normal')

            self.t_pwv.config(state='disable')
            self.e_pwv.config(state='disable')
            #self.t_signalnoise.config(state='disable')
            #self.e_signalnoise.config(state='disable')
            self.t_waverange.config(state='disable')
            self.t_magrange.config(state='normal')
            self.t_kmag.config(state='disable')
            self.e_kmag.config(state='disable')
            self.t_temperature.config(state='disable')
            self.e_temperature.config(state='disable')
            self.wavemin.ShowEntry("disable")
            self.wavemax.ShowEntry("disable")    	    
            
            self.mag.set(" ")
            self.plot_H.config(state='disable')
            self.plot_K.config(state='disable')

            self.t_emission.config(state='disable')
            self.t_restwave.config(state='disable')
            self.e_restwave.config(state='disable')
            #self.u_restwave.config(state='disable')
            self.t_doppler.config(state='disable')
            self.e_doppler.config(state='disable')
            #self.u_doppler.config(state='disable')
            self.t_lineflux.config(state='disable')
            self.e_lineflux.config(state='disable')
            #self.u_lineflux.config(state='disable')
            self.t_linewidth.config(state='disable')
            self.e_linewidth.config(state='disable')
            #self.u_linewidth.config(state='disable')
            
               
        #elif self.mode.get() == "magnitude":  
            #for t, b in self.MAGNITUDE :
            #    t.config(state='active')
            #    b.ShowEntry("active")
                        
            #for t, b in self.SINGLE :
             #   t.config(state='disable')
             #   b.ShowEntry("disable")
            
            #for f in self.PLOT :
            #    f.config(state='disable')
            
            #self.t_pwv.config(state='disable')
            #self.e_pwv.config(state='disable')
            #self.t_signalnoise.config(state='normal')
            #self.e_signalnoise.config(state='normal')
            #self.t_waverange.config(state='disable')
            #self.t_magrange.config(state='disable')
            #self.t_kmag.config(state='disable')
            #self.e_kmag.config(state='disable')
            #self.t_temperature.config(state='disable')
            #self.e_temperature.config(state='disable')
            #self.wavemin.ShowEntry("disable")
            #self.wavemax.ShowEntry("disable")
                                    
            #self.mag.set(" ")
            #self.plot_H.config(state='disable')
            #self.plot_K.config(state='disable')
              
		
        elif self.mode.get() == "plot_wave":
    	    #for t, b in self.MAGNITUDE :
            #    t.config(state='disable')
            #    b.ShowEntry("disable")
    			
            for t in self.SINGLE :
                t.config(state='normal')

    	    #for t, b in self.SINGLE :
            #    t.config(state="disable")
            #    b.ShowEntry("disable")
                
            for f in self.PLOT :
                f.config(state='disable')
    		
            self.t_pwv.config(state='active')
            self.e_pwv.config(state='normal')
            #self.t_signalnoise.config(state='disable')
            #self.e_signalnoise.config(state='disable')
            self.t_waverange.config(state='normal')
            self.t_magrange.config(state='disable')
            self.t_kmag.config(state='normal')
            self.e_kmag.config(state='normal')
            self.t_temperature.config(state='normal')
            self.e_temperature.config(state='normal')

            #self.t_point.config(state='active')
       	    #self.magband.ShowEntry("active")

            self.wavemin.ShowEntry("active")
            self.wavemax.ShowEntry("active")
            
            self.mag.set("plot_wave_H")
            self.MagSelect()
            self.plot_H.config(state='active')
            self.plot_K.config(state='active')

            self.t_emission.config(state='normal')
            self.t_restwave.config(state='normal')
            self.e_restwave.config(state='normal')
            #self.u_restwave.config(state='normal')
            self.t_doppler.config(state='normal')
            self.e_doppler.config(state='normal')
            #self.u_doppler.config(state='normal')
            self.t_lineflux.config(state='normal')
            self.e_lineflux.config(state='normal')
            #self.u_lineflux.config(state='normal')
            self.t_linewidth.config(state='normal')
            self.e_linewidth.config(state='normal')
            #self.u_linewidth.config(state='normal')

        elif self.mode.get() == "plot_signal":
    	    #for t, b in self.MAGNITUDE :
            #    t.config(state='disable')
            #    b.ShowEntry("disable")
    			
            for t in self.SINGLE :
                t.config(state='normal')

    	    #for t, b in self.SINGLE :
            #    t.config(state="disable")
            #    b.ShowEntry("disable")
                
            for f in self.PLOT :
                f.config(state='disable')
    		
            self.t_pwv.config(state='active')
            self.e_pwv.config(state='normal')
            #self.t_signalnoise.config(state='disable')
            #self.e_signalnoise.config(state='disable')
            self.t_waverange.config(state='normal')
            self.t_magrange.config(state='disable')
            self.t_kmag.config(state='normal')
            self.e_kmag.config(state='normal')
            self.t_temperature.config(state='normal')
            self.e_temperature.config(state='normal')

            #self.t_point.config(state='active')
       	    #self.magband.ShowEntry("active")

            self.wavemin.ShowEntry("active")
            self.wavemax.ShowEntry("active")
            
            self.mag.set("plot_wave_H")
            self.MagSelect()
            self.plot_H.config(state='active')
            self.plot_K.config(state='active')

            self.t_emission.config(state='normal')
            self.t_restwave.config(state='normal')
            self.e_restwave.config(state='normal')
            #self.u_restwave.config(state='normal')
            self.t_doppler.config(state='normal')
            self.e_doppler.config(state='normal')
            #self.u_doppler.config(state='normal')
            self.t_lineflux.config(state='normal')
            self.e_lineflux.config(state='normal')
            #self.u_lineflux.config(state='normal')
            self.t_linewidth.config(state='normal')
            self.e_linewidth.config(state='normal')
            #self.u_linewidth.config(state='normal')
            
        elif self.mode.get() == "plot_noise":
    	    #for t, b in self.MAGNITUDE :
            #    t.config(state='disable')
            #    b.ShowEntry("disable")
    			
            for t in self.SINGLE :
                t.config(state='normal')

    	    #for t, b in self.SINGLE :
            #    t.config(state="disable")
            #    b.ShowEntry("disable")
                
            for f in self.PLOT :
                f.config(state='disable')
    		
            self.t_pwv.config(state='active')
            self.e_pwv.config(state='normal')
            #self.t_signalnoise.config(state='disable')
            #self.e_signalnoise.config(state='disable')
            self.t_waverange.config(state='normal')
            self.t_magrange.config(state='disable')
            self.t_kmag.config(state='normal')
            self.e_kmag.config(state='normal')
            self.t_temperature.config(state='normal')
            self.e_temperature.config(state='normal')

            #self.t_point.config(state='active')
       	    #self.magband.ShowEntry("active")

            self.wavemin.ShowEntry("active")
            self.wavemax.ShowEntry("active")
            
            self.mag.set("plot_wave_H")
            self.MagSelect()
            self.plot_H.config(state='active')
            self.plot_K.config(state='active')

            self.t_emission.config(state='normal')
            self.t_restwave.config(state='normal')
            self.e_restwave.config(state='normal')
            #self.u_restwave.config(state='normal')
            self.t_doppler.config(state='normal')
            self.e_doppler.config(state='normal')
            #self.u_doppler.config(state='normal')
            self.t_lineflux.config(state='normal')
            self.e_lineflux.config(state='normal')
            #self.u_lineflux.config(state='normal')
            self.t_linewidth.config(state='normal')
            self.e_linewidth.config(state='normal')
            #self.u_linewidth.config(state='normal')            
            
    def MagSelect(self):
        if self.mag.get() == "plot_wave_H":
            self.wavemin.ActiveEntry('h')
            self.wavemax.ActiveEntry('h')

        elif self.mag.get() == "plot_wave_K":
            self.wavemin.ActiveEntry('k')
            self.wavemax.ActiveEntry('k')
    			
    	    
    def Get_Mag(self):
        kmag = self.e_kmag.get()
        temperature = self.e_temperature.get()
        v_kmag, v_temperature = float(kmag), float(temperature)
        v_mag = Cal_Template_Mag_band(v_kmag, v_temperature)
	#v_mag = ['%.2f' % v_mag[0], '%.2f' % v_mag[1]]
	#self.magband.ShowValue(v_mag)
        self.x.set('%.2f' % (v_mag[0],))
        self.y.set('%.2f' % (v_mag[1],))
	
        

    def Run(self):
        #RESOLUTION = [v_res_j, v_res_h, v_res_k, v_res_l, v_res_m] = self.resband.GetEntry()
        exptime = self.exptime.get()
        expnumber = self.expnumber.get()
        #signoise = self.e_signalnoise.get()
        kmag = self.e_kmag.get()
        temperature = self.e_temperature.get()
        v_exptime, v_expnumber, v_kmag, v_temperature = float(exptime), float(expnumber), float(kmag), float(temperature)  # v_signoise = float(signoise)
        restwave_type = self.e_restwave.get()
        flux_type = self.e_lineflux.get()
        linewidth_type = self.e_linewidth.get()
        doppler_type = self.e_doppler.get()
        v_restwave  = float(restwave_type)   
        v_lineflux  = float(flux_type) * 1E-18
        v_doppler = float(doppler_type)

        # Print ouput
        o_lineflux, o_linewidth = float(flux_type), float(linewidth_type)

        # Change unit to Km s-1
        v_linewidth = (v_restwave * float(linewidth_type)) / (299792500 * pow(10, -3)) # Unit in um
        o_doppler = (v_restwave * float(doppler_type)) / (299792500 * pow(10, -3)) # Unit in um

        # Derive for calibrate line signal peak and line equivalent width based on input of line width and line flux
        v_peaksignal = (v_lineflux * 2.35482) / (np.sqrt(2.0*np.pi)* v_linewidth)
        v_ewidth = (np.sqrt(2.0*np.pi)* v_linewidth) / 2.35482
            
        if self.mode.get() == "single":
            
            v_snr = Mode_SN_band(v_exptime, v_expnumber, v_kmag, v_temperature, v_lineflux, v_linewidth, v_restwave, v_doppler)
            v_mag = Cal_Template_Mag_band(v_kmag, v_temperature)
            
            #self.snrband.ShowValue(v_snr)

            print(' ')
            print('==========================================================================')
            print(ETC_version + ' - The calculation Signal-to-Noise from single magnitude input')
            print(' ')
            print('Exposure Time [s] = %d' %v_exptime)
            print('Exposure Number = %d' %v_expnumber)
            print('K Magnitude [mag] = %d' %v_kmag)
            print('Blackbody Temperature [K] = %d' %v_temperature)
            print('Magnitude values:')
            print('H [mag] = %.2f' %v_mag[0])
            print('K [mag] = %.2f' %v_mag[1])
            print('Signal-to-Noise values:')
            print('H = ', v_snr[0])
            print('K = ', v_snr[1])
            print(' ')
            print('===========================================================================')
            print(' ')
           
        elif self.mode.get() == "plot":
            mag_min = self.e_mag_min.get()
            mag_max = self.e_mag_max.get()
            v_mag_min, v_mag_max = float(mag_min), float(mag_max)
            
            print(' ')
            print('================================================================================')
            print(ETC_version + ' - The calculation Signal-to-Noise vs. Magnitude')
            print(' ')
            print('Exposure Time [s] = %d' %v_exptime)
            print('Exposure Number = %d' %v_expnumber)
            print('Magnitude Range [mag] : %d' %v_mag_min + ' - ' + '%d' %v_mag_max)
            print(' ')
            print('================================================================================')
            print(' ')

            Mode_SN_Mag_Plot(v_exptime, v_expnumber, v_mag_max, v_mag_min, v_lineflux, v_linewidth, v_restwave, v_doppler)        
            
        elif self.mode.get() == "magnitude":
            v_magnitude = Mode_Mag_Limiting_band(v_exptime, v_expnumber, v_signoise)  
            
            self.signoiseband.ShowValue(v_magnitude)  

            print(' ')
            print('===============================================================================')
            print(ETC_version + ' - The calculation magnitude values from signal-to-noise input')
            print(' ')
            print('Exposure Time [s] = %d' %v_exptime)
            print('Exposure Number = %d' %v_expnumber)
            print('Signal-to-Noise = %d' %v_signoise)
            print('Magnitude values:')
            print('H [mag] = ', v_magnitude[0])
            print('K [mag] = ', v_magnitude[1])
            print(' ')
            print('===============================================================================')
            print(' ')

        elif self.mode.get() == "plot_wave":

            if self.mag.get() == "plot_wave_H":            
                #MAGNITUDE = [v_mag_j, v_mag_h, v_mag_k, v_mag_l, v_mag_m] = self.magband.GetEntry()      
                pwv_type = self.e_pwv.get()   
                waveminH = self.wavemin.GetActiveEntryValue("h")
                wavemaxH = self.wavemin.GetActiveEntryValue("k") #self.wavemax.GetActiveEntryValue("h")            
                v_pwv, v_waveminH, v_wavemaxH = float(pwv_type), float(waveminH), float(wavemaxH)		

                print(' ')
                print('==========================================================================')
                print(ETC_version + ' - The calculation Signal-to-Noise vs. Wavelength')
                print(' ')
                print('Exposure Time [s] = %d' % v_exptime)
                print('Exposure Number = %d' % v_expnumber)
                print('K Magnitude [mag] = %d' %v_kmag)
                print('Blackbody Temperature [K] = %d' %v_temperature)
                print('PWV = %d' %v_pwv)
                print('Calculated Wavelength Range [um] : ' + '%0.6f' %v_waveminH + ' - ' + '%0.6f' %v_wavemaxH)
                print('Rest Wavelength [um] = %0.6f' %v_restwave)
                print('Line Flux [1E-18 W m-2] = %d' %o_lineflux)
                print('Line Width [um] = %0.6f' %v_linewidth)
                print('Calibrate Signal Line Peak [W m-2 um-1] = %0.3e' %v_peaksignal)
                print('Line Equivalent Width [um] = %0.6f' %v_ewidth)
                print('Doppler Shift [um] = %0.6f' %o_doppler)
                print(' ')
                print('==========================================================================')
                print(' ')
                
                Mode_SN_wave_Plot(v_exptime, v_expnumber, v_pwv, v_waveminH, v_wavemaxH, v_kmag, v_temperature, v_lineflux \
                                   , v_linewidth, v_restwave, v_doppler)

            elif self.mag.get() == "plot_wave_K":            
                #MAGNITUDE = [v_mag_j, v_mag_h, v_mag_k, v_mag_l, v_mag_m] = self.magband.GetEntry()      
                pwv_type = self.e_pwv.get()   
                waveminK = self.wavemax.GetActiveEntryValue("h") #self.wavemin.GetActiveEntryValue("k")
                wavemaxK = self.wavemax.GetActiveEntryValue("k")           
                v_pwv, v_waveminK, v_wavemaxK = float(pwv_type), float(waveminK), float(wavemaxK)			

                print(' ')
                print('==========================================================================')
                print(ETC_version + ' - The calculation Signal-to-Noise vs. Wavelength')
                print(' ')
                print('Exposure Time [s] = %d' %v_exptime)
                print('Exposure Number = %d' %v_expnumber)
                print('K Magnitude [mag] = %d' %v_kmag)
                print('Blackbody Temperature [K] = %d' %v_temperature)
                print('PWV = %d' %v_pwv)
                print('Calculated Wavelength Range [um] : %0.6f' %v_waveminK + ' - ' + '%0.6f' %v_wavemaxK)
                print('Rest Wavelength [um] = %0.6f' %v_restwave)
                print('Line Flux [1E-18 W m-2] = %d' %o_lineflux)
                print('Line Width [um] = %0.6f' %v_linewidth)
                print('Calibrate Signal Line Peak [W m-2 um-1] = %0.3e' %v_peaksignal)
                print('Line Equivalent Width [um] = %0.6f' %v_ewidth)
                print('Doppler Shift [um] = %0.6f' % o_doppler)
                print(' ')
                print('==========================================================================')
                print(' ')

                Mode_SN_wave_Plot(v_exptime, v_expnumber, v_pwv, v_waveminK, v_wavemaxK, v_kmag, v_temperature, v_lineflux \
                                   , v_linewidth, v_restwave, v_doppler)

        elif self.mode.get() == "plot_signal":
            
            if self.mag.get() == "plot_wave_H":            
                #MAGNITUDE = [v_mag_j, v_mag_h, v_mag_k, v_mag_l, v_mag_m] = self.magband.GetEntry()      
                pwv_type = self.e_pwv.get()   
                waveminH = self.wavemin.GetActiveEntryValue("h")
                wavemaxH = self.wavemin.GetActiveEntryValue("k") #self.wavemax.GetActiveEntryValue("h")            
                v_pwv, v_waveminH, v_wavemaxH = float(pwv_type), float(waveminH), float(wavemaxH)
                
                print(' ')
                print('==========================================================================')
                print(ETC_version + ' - The calculation Signal with Noise vs. Wavelength')
                print(' ')
                print('Exposure Time [s] = %d' %v_exptime)
                print('Exposure Number = %d' %v_expnumber)
                print('K Magnitude [mag] = %d' %v_kmag)
                print('Blackbody Temperature [K] = %d' %v_temperature)
                print('PWV = %d' %v_pwv)
                print('Calculated Wavelength Range [um] : %0.6f' %v_waveminH + ' - ' + '%0.6f' % v_wavemaxH)
                print('Rest Wavelength [um] = %0.6f' %v_restwave)
                print('Line Flux [1E-18 W m-2] = %d' %o_lineflux)
                print('Line Width [um] = %0.6f' %v_linewidth)
                print('Calibrate Signal Line Peak [W m-2 um-1] = %0.3e' %v_peaksignal)
                print('Line Equivalent Width [um] = %0.6f' % v_ewidth)
                print('Doppler Shift [um] = %0.6f' % o_doppler)
                
                #print '=========================================================================='

                #Mode_MS_wave_Plot(v_exptime, v_expnumber, v_pwv, v_waveminH, v_wavemaxH, v_kmag, v_temperature)
                Mode_MSL_wave_Plot(v_exptime, v_expnumber, v_pwv, v_waveminH, v_wavemaxH, v_kmag, v_temperature, v_lineflux \
                                   , v_linewidth, v_restwave, v_doppler)

            elif self.mag.get() == "plot_wave_K":            
                #MAGNITUDE = [v_mag_j, v_mag_h, v_mag_k, v_mag_l, v_mag_m] = self.magband.GetEntry()      
                pwv_type = self.e_pwv.get()   
                waveminK = self.wavemax.GetActiveEntryValue("h") #self.wavemin.GetActiveEntryValue("k")
                wavemaxK = self.wavemax.GetActiveEntryValue("k")                   
                v_pwv, v_waveminK, v_wavemaxK = float(pwv_type), float(waveminK), float(wavemaxK)			

                print(' ')
                print('==========================================================================')
                print(ETC_version + ' - The calculation Signal with Noise vs. Wavelength')
                print(' ')
                print('Exposure Time [s] = %d' %v_exptime)
                print('Exposure Number = %d' %v_expnumber)
                print('K Magnitude [mag] = %d' %v_kmag)
                print('Blackbody Temperature [K] = %d' %v_temperature)
                print('PWV = %d' %v_pwv)
                print('Calculated Wavelength Range [um] : %0.6f' %v_waveminK + ' - ' + '%0.6f' %v_wavemaxK)
                print('Rest Wavelength [um] = %0.6f' %v_restwave)
                print('Line Flux [1E-18 W m-2] = %d' %o_lineflux)
                print('Line Width [um] = %0.6f' %v_linewidth)
                print('Calibrate Signal Line Peak [W m-2 um-1] = %0.3e' %v_peaksignal)
                print('Line Equivalent Width [um] = %0.6f' %v_ewidth)
                print('Doppler Shift [um] = %0.6f' %o_doppler)
                
                #print '=========================================================================='

                #Mode_MS_wave_Plot(v_exptime, v_expnumber, v_pwv, v_waveminK, v_wavemaxK, v_kmag, v_temperature)    
                Mode_MSL_wave_Plot(v_exptime, v_expnumber, v_pwv, v_waveminK, v_wavemaxK, v_kmag, v_temperature, v_lineflux \
                                   , v_linewidth, v_restwave, v_doppler)


        elif self.mode.get() == "plot_noise":

            if self.mag.get() == "plot_wave_H":            
                #MAGNITUDE = [v_mag_j, v_mag_h, v_mag_k, v_mag_l, v_mag_m] = self.magband.GetEntry()      
                pwv_type = self.e_pwv.get()   
                waveminH = self.wavemin.GetActiveEntryValue("h")
                wavemaxH = self.wavemin.GetActiveEntryValue("k") #self.wavemax.GetActiveEntryValue("h")          
                v_pwv, v_waveminH, v_wavemaxH = float(pwv_type), float(waveminH), float(wavemaxH)
                
                print(' ')
                print('==========================================================================')
                print(ETC_version + ' - The calculation Noise vs. Wavelength')
                print(' ')
                print('Exposure Time [s] = %d' %v_exptime)
                print('Exposure Number = %d' %v_expnumber)
                print('K Magnitude [mag] = %d' %v_kmag)
                print('Blackbody Temperature [K] = %d' %v_temperature)
                print('PWV = %d' %v_pwv)
                print('Calculated Wavelength Range [um] : %0.6f' %v_waveminH + ' - ' + '%0.6f' %v_wavemaxH)
                print('Rest Wavelength [um] = %0.6f' %v_restwave)
                print('Line Flux [1E-18 W m-2] = %d' %o_lineflux)
                print('Line Width [um] = %0.6f' %v_linewidth)
                print('Calibrate Signal Line Peak [W m-2 um-1] = %0.3e' %v_peaksignal)
                print('Line Equivalent Width [um] = %0.6f' %v_ewidth)
                print('Doppler Shift [um] = %0.6f' %o_doppler)
                print(' ')
                print('==========================================================================')

                #Mode_MN_wave_Plot(v_exptime, v_expnumber, v_pwv, v_waveminH, v_wavemaxH, v_kmag, v_temperature)
                Mode_MN_wave_Plot(v_exptime, v_expnumber, v_pwv, v_waveminH, v_wavemaxH, v_kmag, v_temperature, v_lineflux \
                                   , v_linewidth, v_restwave, v_doppler)

            elif self.mag.get() == "plot_wave_K":            
                #MAGNITUDE = [v_mag_j, v_mag_h, v_mag_k, v_mag_l, v_mag_m] = self.magband.GetEntry()      
                pwv_type = self.e_pwv.get()   
                waveminK = self.wavemax.GetActiveEntryValue("h") #self.wavemin.GetActiveEntryValue("k")
                wavemaxK = self.wavemax.GetActiveEntryValue("k")                 
                v_pwv, v_waveminK, v_wavemaxK = float(pwv_type), float(waveminK), float(wavemaxK)			

                print(' ')
                print('==========================================================================')
                print(ETC_version + ' - The calculation Noise vs. Wavelength')
                print(' ')
                print('Exposure Time [s] = %d' %v_exptime)
                print('Exposure Number = %d' %v_expnumber)
                print('K Magnitude [mag] = %d' %v_kmag)
                print('Blackbody Temperature [K] = %d' %v_temperature)
                print('PWV = %d' %v_pwv)
                print('Calculated Wavelength Range [um] : %0.6f' %v_waveminK + ' - ' + '%0.6f' % v_wavemaxK)
                print('Rest Wavelength [um] = %0.6f' %v_restwave)
                print('Line Flux [1E-18 W m-2] = %d' %o_lineflux)
                print('Line Width [um] = %0.6f' %v_linewidth)
                print('Calibrate Signal Line Peak [W m-2 um-1] = %0.3e' %v_peaksignal)
                print('Line Equivalent Width [um] = %0.6f' %v_ewidth)
                print('Doppler Shift [um] = %0.6f' %o_doppler)
                print(' ')
                print('==========================================================================')
                print(' ')

                #Mode_MN_wave_Plot(v_exptime, v_expnumber, v_pwv, v_waveminK, v_wavemaxK, v_kmag, v_temperature)    
                Mode_MN_wave_Plot(v_exptime, v_expnumber, v_pwv, v_waveminK, v_wavemaxK, v_kmag, v_temperature, v_lineflux \
                                   , v_linewidth, v_restwave, v_doppler)        

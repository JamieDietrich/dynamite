### DYNAmical Multi-planet Injection TEster (DYNAMITE) ###
### GUI ###
### Jeremy Dietrich ###
### jdietrich1@email.arizona.edu ###
### 2022 January 7 ###
### Version 2.0 ###
### Dietrich & Apai (2020), AJ, 160, 107D ###
### Dietrich & Apai (2021), AJ, 161, 17D ###
### Dietrich, Apai, & Malhotra (2022), accepted to AJ ###
### https://iopscience.iop.org/article/10.3847/1538-3881/aba61d ###

import os
import sys
import ast
import time
import atexit
import socket
import threading
import subprocess
import numpy as np
from tkinter import ttk
import tkinter as Tkinter
import dynamite_targets_editor
from tkinter import filedialog as tkFD
from tkinter import messagebox as tkMB
import tkinter.scrolledtext as ScrolledText
from tkinter import constants as Tkconstants
from dynamite_targets_v2 import dynamite_targets

class CreateToolTip(object):
    """Modified from https://stackoverflow.com/questions/3221956/how-do-i-display-tooltips-in-tkinter"""
    def __init__(self, widget, root, tooltips, config_entries, widget1 = None):
        self.root = root
        self.widget = widget
        self.tw = None
        if config_entries != None:
            try:   self.text = tooltips[[k for k,v in config_entries.items() if v == (widget if widget1 == None else widget1)][0]]
            except:self.text = "No tooltip entry found in config file"
        else:
            self.text = tooltips
        self.widget.bind("<ButtonPress-3>", self.enter)
        self.widget.bind("<Leave>", self.hidetip)
        self.widget.bind("<ButtonPress>", self.hidetip)
        
    def enter(self, event = None):
        if self.tw == None:
            self.widget.after(7000 if len(self.text) > 150 else 7000 if len(self.text) > 100 else 5000, self.hidetip)
            self.showtip()

    def showtip(self, event = None):
        x, y, _cx, _cy = self.widget.bbox("insert")
        x += self.widget.winfo_rootx() + 25
        y += self.widget.winfo_rooty() + 30
        self.tw = Tkinter.Toplevel(self.root)
        self.tw.wm_overrideredirect(True)
        self.tw.wm_geometry("+%d+%d" % (x, y))
        label = ttk.Label(self.tw, text=self.text, justify='left',background="#ffffff", relief='solid', borderwidth=1,wraplength = 600)
        label.pack(ipadx=1)

    def hidetip(self, event = None):
        if self.tw:
            self.tw.destroy()
            self.tw = None

class DynamiteGUI:
    def __init__(self):
        """Creates an instance of the GUI with instance variables."""                   
        self.config_entries = {}
        self.run_editors = {}
        self.plot_editors = {}
        self.root = Tkinter.Tk()
        self.root.minsize(1200, 400)
                
        if (len(sys.argv)  >=  2):
            self.load_config_data(sys.argv[1])
        else:
            self.load_config_data("dynamite_config.txt")
            
        ttk.Style().map("TCombobox", fieldbackground = [('!readonly', '#FFFFFF'),('readonly', '#FFFFFF'),],)
        self.root.title("DYNAMITE ---> " + self.config_file + " (right-click on field for tooltips)")
                        
        self.mc_chain_values = "10000"
        self.system_values = sorted(list(dynamite_targets(self.config_parameters).targets.keys()),key=str.casefold)
        self.radmin_values = "0.1"
        self.radmax_values = "5"
        self.plot_colors_values = '["377eb8", "d95f02", "4daf4a", "984ea3", "e41a1c", "e6ab02"]'
        self.additional_values = "[[]]"
        self.unconfirmed_values = "[[]]"
        self.removed_values = "[]" 
        self.mode_values = ["single", "tess", "kepler", "k2", "all", "test"]
        self.period_values = ["epos", "syssim"]
        self.radius_values = ["epos", "syssim"]
        self.radtype_values = ["powerlaw", "clustered"]
        self.mass_radius_values = ["mrexo", "otegi"]
        self.otegi_rho_values = ["rocky", "volatile"]
        self.inclination_values =["rayleigh_iso", "syssim"]
        self.eccentricity_values =["rayleigh", "syssim"]
        self.stability_values = ["hill", "spock", "specfrac"]
        self.plot_values = [True, False]
        self.show_plots_values = [True, False]
        self.save_plots_values = [True, False]
        self.plot_p_r_values = [True, False]
        self.plot_p_scale_values = ["linear", "log"]
        self.plot_tdtp_values = [True, False]
        self.plot_deltas_values = [True, False]
        self.plot_ratios_values = [True, False]
        self.plot_hist_values = [True, False]
        self.plot_pdf_values = [True, False]
        self.show_rearth_values = [True, False]
        self.use_mass_values = [True, False]
        self.add_unconfirmed_values = [True, False]
        self.ind_p_values = ["linear_zoom", "linear", "log"]
        self.ind_r_values = ["linear_zoom", "linear"]
        self.ind_i_values = ["full", "truncated"]
        
        self.mc_chain = Tkinter.StringVar()
        self.mc_chain.set(self.mc_chain_values)
        self.config_entries["MC_chain"] = self.mc_chain
        
        self.radmin = Tkinter.StringVar()
        self.radmin.set(self.radmin_values)
        self.config_entries["radmin"] = self.radmin
        
        self.radmax = Tkinter.StringVar()
        self.radmax.set(self.radmax_values)
        self.config_entries["radmax"] = self.radmax
        
        self.additional = Tkinter.StringVar()
        self.additional.set(self.additional_values)
        self.config_entries["additional"] = self.additional
        
        self.unconfirmed = Tkinter.StringVar()
        self.unconfirmed.set(self.unconfirmed_values)
        self.config_entries["unconfirmed"] = self.unconfirmed
        
        self.removed = Tkinter.StringVar()
        self.removed.set(self.removed_values)
        self.config_entries["removed"] = self.removed
        
        self.plot_colors = Tkinter.StringVar()
        self.plot_colors.set(self.plot_colors_values)
        self.config_entries["plt_colors"] = self.plot_colors

        self.mode_box = ttk.Combobox(self.root, width = 10, justify = Tkconstants.RIGHT)
        self.mode_box["values"] = self.mode_values
        self.mode_box.current(0)
        self.config_entries["mode"] = self.mode_box
        self.mode_box.bind('<<ComboboxSelected>>', self.set_plot_save_false)
        self.mode_box.configure(state = 'readonly')
        
        self.system_box = ttk.Combobox(self.root, width = 18, height = 20, justify = Tkconstants.RIGHT)
        self.system_box["values"] = self.system_values
        self.system_box.current(0)
        self.config_entries["system"] = self.system_box
        self.system_box.bind('<<ComboboxSelected>>', self.set_plot_save_false)
        self.system_box.configure(state = 'readonly')
        
        self.period_box = ttk.Combobox(self.root, width = 10, justify = Tkconstants.RIGHT)
        self.period_box["values"] = self.period_values
        self.period_box.current(0)
        self.config_entries["period"] = self.period_box
        self.period_box.bind('<<ComboboxSelected>>', self.set_plot_save_false)
        self.period_box.configure(state = 'readonly')
        
        self.radius_box = ttk.Combobox(self.root, width = 10, justify = Tkconstants.RIGHT)
        self.radius_box["values"] = self.radius_values
        self.radius_box.current(1)
        self.config_entries["radius"] = self.radius_box
        self.radius_box.bind('<<ComboboxSelected>>', self.set_plot_save_false)
        self.radius_box.configure(state = 'readonly')
        
        self.radtype_box = ttk.Combobox(self.root, width = 11, justify = Tkconstants.RIGHT)
        self.radtype_box["values"] = self.radtype_values
        self.radtype_box.current(1)
        self.config_entries["radtype"] = self.radtype_box
        self.radtype_box.bind('<<ComboboxSelected>>', self.set_plot_save_false)
        self.radtype_box.configure(state = 'readonly')

        self.mass_radius_box = ttk.Combobox(self.root, width = 10, justify = Tkconstants.RIGHT)
        self.mass_radius_box["values"] = self.mass_radius_values
        self.mass_radius_box.current(1)
        self.config_entries["mass_radius"] = self.mass_radius_box
        self.mass_radius_box.bind('<<ComboboxSelected>>', self.set_plot_save_false)
        self.mass_radius_box.configure(state = 'readonly')
        
        self.otegi_rho_box = ttk.Combobox(self.root, width = 10, justify = Tkconstants.RIGHT)
        self.otegi_rho_box["values"] = self.otegi_rho_values
        self.otegi_rho_box.current(1)
        self.config_entries["otegi_rho"] = self.otegi_rho_box
        self.otegi_rho_box.bind('<<ComboboxSelected>>', self.set_plot_save_false)
        self.otegi_rho_box.configure(state = 'readonly')
        
        self.inclination_box = ttk.Combobox(self.root, width = 12, justify = Tkconstants.RIGHT)
        self.inclination_box["values"] = self.inclination_values
        self.inclination_box.current(1)
        self.config_entries["inclination"] = self.inclination_box
        self.inclination_box.bind('<<ComboboxSelected>>', self.set_plot_save_false)
        self.inclination_box.configure(state = 'readonly')
        
        self.eccentricity_box = ttk.Combobox(self.root, width = 8, justify = Tkconstants.RIGHT)
        self.eccentricity_box["values"] = self.eccentricity_values
        self.eccentricity_box.current(1)
        self.config_entries["eccentricity"] = self.eccentricity_box
        self.eccentricity_box.bind('<<ComboboxSelected>>', self.set_plot_save_false)
        self.eccentricity_box.configure(state = 'readonly')
 
        self.stability_box = ttk.Combobox(self.root, width = 10, justify = Tkconstants.RIGHT)
        self.stability_box["values"] = self.stability_values
        self.stability_box.current(1)
        self.config_entries["stability"] = self.stability_box
        self.stability_box.bind('<<ComboboxSelected>>', self.set_plot_save_false)
        self.stability_box.configure(state = 'readonly')
               
        self.plot_box = ttk.Combobox(self.root, width = 10, justify = Tkconstants.RIGHT)
        self.plot_box["values"] = self.plot_values
        self.plot_box.current(1)
        self.config_entries["plot"] = self.plot_box
        self.plot_box.bind('<<ComboboxSelected>>')
        self.plot_box.configure(state = 'readonly')
        
        self.show_plots_box = ttk.Combobox(self.root, width = 10, justify = Tkconstants.RIGHT)
        self.show_plots_box["values"] = self.show_plots_values
        self.show_plots_box.current(1)
        self.config_entries["show_plots"] = self.show_plots_box
        self.show_plots_box.bind('<<ComboboxSelected>>')
        self.show_plots_box.configure(state = 'readonly')
        
        self.save_plots_box = ttk.Combobox(self.root, width = 10, justify = Tkconstants.RIGHT)  
        self.save_plots_box["values"] = self.save_plots_values
        self.save_plots_box.current(1)
        self.config_entries["saved"] = self.save_plots_box
        self.save_plots_box.bind('<<ComboboxSelected>>')
        self.save_plots_box.configure(state = 'readonly')
        
        self.plot_p_r_box = ttk.Combobox(self.root, width = 10, justify = Tkconstants.RIGHT)
        self.plot_p_r_box["values"] = self.plot_p_r_values
        self.plot_p_r_box.current(1)
        self.config_entries["plt_P_R"] = self.plot_p_r_box
        self.plot_p_r_box.bind('<<ComboboxSelected>>')
        self.plot_p_r_box.configure(state = 'readonly')
        
        self.plot_p_scale_box = ttk.Combobox(self.root, width = 10, justify = Tkconstants.RIGHT)
        self.plot_p_scale_box["values"] = self.plot_p_scale_values
        self.plot_p_scale_box.current(1)
        self.config_entries["plt_P_scale"] = self.plot_p_scale_box
        self.plot_p_scale_box.bind('<<ComboboxSelected>>')
        self.plot_p_scale_box.configure(state = 'readonly')
        
        self.plot_tdtp_box = ttk.Combobox(self.root, width = 10, justify = Tkconstants.RIGHT)
        self.plot_tdtp_box["values"] = self.plot_tdtp_values
        self.plot_tdtp_box.current(1)
        self.config_entries["plt_tdtp"] = self.plot_tdtp_box
        self.plot_tdtp_box.bind('<<ComboboxSelected>>')
        self.plot_tdtp_box.configure(state = 'readonly')
        
        self.plot_deltas_box = ttk.Combobox(self.root, width = 10, justify = Tkconstants.RIGHT)
        self.plot_deltas_box["values"] = self.plot_deltas_values
        self.plot_deltas_box.current(1)
        self.config_entries["plt_deltas"] = self.plot_deltas_box
        self.plot_deltas_box.bind('<<ComboboxSelected>>')
        self.plot_deltas_box.configure(state = 'readonly')
        
        self.plot_ratios_box = ttk.Combobox(self.root, width = 10, justify = Tkconstants.RIGHT)
        self.plot_ratios_box["values"] = self.plot_ratios_values
        self.plot_ratios_box.current(1)
        self.config_entries["plt_ratios"] = self.plot_ratios_box
        self.plot_ratios_box.bind('<<ComboboxSelected>>')
        self.plot_ratios_box.configure(state = 'readonly')
        
        self.plot_hist_box = ttk.Combobox(self.root, width = 10, justify = Tkconstants.RIGHT)
        self.plot_hist_box["values"] = self.plot_hist_values
        self.plot_hist_box.current(0)
        self.config_entries["plt_indpars"] = self.plot_hist_box
        self.plot_hist_box.bind('<<ComboboxSelected>>')
        self.plot_hist_box.configure(state = 'readonly')
        
        self.plot_pdf_box = ttk.Combobox(self.root, width = 10, justify = Tkconstants.RIGHT)
        self.plot_pdf_box["values"] = self.plot_pdf_values
        self.plot_pdf_box.current(0)
        self.config_entries["plt_PDFs"] = self.plot_pdf_box
        self.plot_pdf_box.bind('<<ComboboxSelected>>')
        self.plot_pdf_box.configure(state = 'readonly')
        
        self.show_rearth_box = ttk.Combobox(self.root, width = 10, justify = Tkconstants.RIGHT)
        self.show_rearth_box["values"] = self.show_rearth_values
        self.show_rearth_box.current(1)
        self.config_entries["show_Rearth"] = self.show_rearth_box
        self.show_rearth_box.bind('<<ComboboxSelected>>')
        self.show_rearth_box.configure(state = 'readonly')
        
        self.use_mass_box = ttk.Combobox(self.root, width = 10, justify = Tkconstants.RIGHT)
        self.use_mass_box["values"] = self.use_mass_values
        self.use_mass_box.current(1)
        self.config_entries["use_mass"] = self.use_mass_box
        self.use_mass_box.bind('<<ComboboxSelected>>')
        self.use_mass_box.configure(state = 'readonly')
                
        self.ind_p_box = ttk.Combobox(self.root, width = 11, justify = Tkconstants.RIGHT)
        self.ind_p_box["values"] = self.ind_p_values
        self.ind_p_box.current(2)
        self.config_entries["ind_P"] = self.ind_p_box
        self.ind_p_box.bind('<<ComboboxSelected>>')
        self.ind_p_box.configure(state = 'readonly')
        
        self.ind_r_box = ttk.Combobox(self.root, width = 11, justify = Tkconstants.RIGHT)
        self.ind_r_box["values"] = self.ind_r_values
        self.ind_r_box.current(1)
        self.config_entries["ind_R"] = self.ind_r_box
        self.ind_r_box.bind('<<ComboboxSelected>>')
        self.ind_r_box.configure(state = 'readonly')
        
        self.ind_i_box = ttk.Combobox(self.root, width = 11, justify = Tkconstants.RIGHT)
        self.ind_i_box["values"] = self.ind_i_values
        self.ind_i_box.current(0)
        self.config_entries["ind_i"] = self.ind_i_box
        self.ind_i_box.bind('<<ComboboxSelected>>')
        self.ind_i_box.configure(state = 'readonly')
        
        self.add_unconfirmed_box = ttk.Combobox(self.root, width = 6, justify = Tkconstants.RIGHT)
        self.add_unconfirmed_box["values"] = self.add_unconfirmed_values
        self.add_unconfirmed_box.current(0)
        self.config_entries["add_unconfirmed"] = self.add_unconfirmed_box
        self.add_unconfirmed_box.bind('<<ComboboxSelected>>')
        self.add_unconfirmed_box.configure(state = 'readonly')
        
        vcmd = (self.root.register(self.validate_key), '%P', '%S')
        self.check_valid_entries_list = [self.additional, self.unconfirmed, self.removed, self.plot_colors]
        
        row = 0
        self.root.grid_rowconfigure(row,weight=1)
        Tkinter.Label(self.root, text = "V 2", font = ("Helvetica", 7)).grid(row = row, column = 5, sticky = 'EN')
        Tkinter.Label(self.root, text = "PREDICTION PARAMETERS").grid(row = row, column = 1, sticky = 'EW', padx = 0, pady = 5, columnspan = 3)
        row += 1
        self.root.grid_rowconfigure(row,weight=1)
        Tkinter.Label(self.root, text = "MC Chain").grid(row = row, column = 0, padx = 5, sticky = 'W', pady = 5)
        self.mc_chain_text_box = Tkinter.Entry(self.root, textvariable = self.mc_chain, width = 10, validate = "key", validatecommand = vcmd, justify ='center')
        self.mc_chain_text_box.grid(row = row, column = 0, padx = 0, sticky = 'E')
        Tkinter.Label(self.root, text = "Mode").grid(row = row, column = 1, padx = 10, sticky = 'W')
        self.mode_box.grid(row = row, column = 1, padx = 0, sticky = 'E')
        Tkinter.Label(self.root, text = "System").grid(row = row, column = 2, padx = 10, sticky = 'W')
        self.system_box.grid(row = row, column = 2, padx = 0, sticky = 'E')
        Tkinter.Label(self.root, text = "Period").grid(row = row, column = 3, padx = 10, sticky = 'W')
        self.period_box.grid(row = row, column = 3, padx = 20, sticky = 'E')
        Tkinter.Label(self.root, text = "Mass Radius").grid(row = row, column = 4, padx = 0, sticky = 'W')
        self.mass_radius_box.grid(row = row, column = 4, padx = 5, sticky = 'E')
        row += 1
        self.root.grid_rowconfigure(row,weight=1)        
        Tkinter.Label(self.root, text = "Radius").grid(row = row, column = 0, padx = 5, sticky = 'W', pady = 5)
        self.radius_box.grid(row = row, column = 0, padx = 0, sticky = 'E')
        Tkinter.Label(self.root, text = "Rad Type").grid(row = row, column = 1, padx = 10, sticky ='W')
        self.radtype_box.grid(row = row, column = 1, padx = 0, sticky = 'E')
        Tkinter.Label(self.root, text = "Rad Min").grid(row = row, column = 2, sticky = 'W', padx = 10)
        self.radmin_text_box = Tkinter.Entry(self.root, textvariable = self.radmin, width = 7, validate = "key", validatecommand = vcmd, justify ='center')
        self.radmin_text_box.grid(row = row, column = 2, padx = 0, sticky = 'E')
        Tkinter.Label(self.root, text = "Rad Max").grid(row = row, column = 3, sticky = 'W', padx = 10)
        self.radmax_text_box = Tkinter.Entry(self.root, textvariable = self.radmax, width = 7, validate = "key", validatecommand = vcmd, justify ='center')        
        self.radmax_text_box.grid(row = row, column = 3, padx = 20, sticky = 'E')
        Tkinter.Label(self.root, text = "Otegi Rho").grid(row = row, column = 4, padx = 0, sticky ='W')
        self.otegi_rho_box.grid(row = row, column = 4, padx = 5, sticky = 'E')
        row += 1
        self.root.grid_rowconfigure(row,weight=1)
        Tkinter.Label(self.root, text = "Inclination").grid(row = row, column = 0, padx = 5, sticky ='W', pady = 5)
        self.inclination_box.grid(row = row, column = 0, padx = 0, sticky = 'E')
        Tkinter.Label(self.root, text = "Eccentricity").grid(row = row, column = 2, padx = 10, sticky ='W', pady = 5)
        self.eccentricity_box.grid(row = row, column = 2, padx = 0, sticky = 'E')
        Tkinter.Label(self.root, text = "Stability").grid(row = row, column = 4, padx = 0, sticky ='W', pady = 5)
        self.stability_box.grid(row = row, column = 4, padx = 0, sticky = 'E')
        row += 1
        self.root.grid_rowconfigure(row,weight=1)
        Tkinter.Label(self.root, text = "PLOTTING PARAMETERS").grid(row = row, column = 1, sticky = 'EW', padx = 0, pady = 5, columnspan = 3)
        row += 1
        self.root.grid_rowconfigure(row,weight=1)
        Tkinter.Label(self.root, text = "Create Plots").grid(row = row, column = 0, padx = 5, sticky = 'W', pady = 5)
        self.plot_box.grid(row = row, column = 0, padx = 10, sticky = 'E')
        Tkinter.Label(self.root, text = "Show Plots").grid(row = row, column = 1, padx = 10, sticky = 'W')
        self.show_plots_box.grid(row = row, column = 1, padx = 0, sticky = 'E')
        Tkinter.Label(self.root, text = "Plot P vs R").grid(row = row, column = 2, sticky = 'W', padx = 10)
        self.plot_p_r_box.grid(row = row, column = 2, padx = 0, sticky = 'E')
        Tkinter.Label(self.root, text = "Use Saved Data").grid(row = row, column = 3, padx = 10, sticky = 'W')
        self.save_plots_box.grid(row = row, column = 3, padx = 0, sticky = 'E')
        Tkinter.Label(self.root, text = "IND P").grid(row = row, column = 4, sticky = 'W', padx = 5)
        self.ind_p_box.grid(row = row, column = 4, padx = 5, sticky = 'E')
        row += 1
        self.root.grid_rowconfigure(row,weight=1)
        Tkinter.Label(self.root, text = "Plot P Scale").grid(row = row, column = 0, padx = 5, sticky = 'W', pady = 5)
        self.plot_p_scale_box.grid(row = row, column = 0, padx = 10, sticky = 'E')
        Tkinter.Label(self.root, text = "Plot TD-TP").grid(row = row, column = 1, padx = 10, sticky = 'W')
        self.plot_tdtp_box.grid(row = row, column = 1, padx = 0, sticky = 'E')
        Tkinter.Label(self.root, text = "Plot Deltas").grid(row = row, column = 2, padx = 10, sticky = 'W')
        self.plot_deltas_box.grid(row = row, column = 2, padx =  0, sticky = 'E')
        Tkinter.Label(self.root, text = "Plot Ratios").grid(row = row, column = 3, sticky = 'W', padx = 10)
        self.plot_ratios_box.grid(row = row, column = 3, padx = 0, sticky = 'E')
        Tkinter.Label(self.root, text = "IND R").grid(row = row, column = 4, padx = 5, sticky = 'W')
        self.ind_r_box.grid(row = row, column = 4, padx = 5, sticky = 'E')
        row += 1
        self.root.grid_rowconfigure(row,weight=1)
        Tkinter.Label(self.root, text = "Plot Hist").grid(row = row, column = 0, padx = 5, sticky = 'W', pady =5)
        self.plot_hist_box.grid(row = row, column = 0, padx = 10, sticky = 'E')
        Tkinter.Label(self.root, text = "Plot PDF").grid(row = row, column = 1, padx = 10, sticky = 'W')
        self.plot_pdf_box.grid(row = row, column = 1, padx = 0, sticky = 'E')
        Tkinter.Label(self.root, text = "Show Rearth").grid(row = row, column = 2, padx = 10, sticky = 'W')
        self.show_rearth_box.grid(row = row, column = 2, padx =0, sticky = 'E')
        Tkinter.Label(self.root, text = "Use Mass").grid(row = row, column = 3, padx = 10, sticky = 'W')
        self.use_mass_box.grid(row = row, column = 3, padx =0, sticky = 'E')
        Tkinter.Label(self.root, text = "IND I").grid(row = row, column = 4, padx = 5, sticky = 'W')
        self.ind_i_box.grid(row = row, column = 4, padx = 5, sticky = 'E')
        row += 1
        self.root.grid_rowconfigure(row,weight=1)
        Tkinter.Label(self.root, text = "Plot Colors").grid(row = row, column = 0, sticky = 'E', padx = 15, pady = 5)
        self.plot_colors_text_box = Tkinter.Entry(self.root, textvariable = self.plot_colors, width = 90, justify ='left')
        self.plot_colors_text_box.grid(row = row, column = 1, padx = 5, sticky = 'EW',columnspan = 4)
        row += 1
        self.root.grid_rowconfigure(row,weight=1)
        Tkinter.Label(self.root, text = "SYSTEM ARCHITECTURE PARAMETERS").grid(row = row, column = 1, sticky = 'EW', padx = 120, pady = 5, columnspan = 3)
        row += 1
        self.root.grid_rowconfigure(row,weight=1)
        Tkinter.Label(self.root, text = "Additional Planets").grid(row = row, column = 0, sticky = 'E', padx = 15, pady = 5)
        self.additional_text_box = Tkinter.Entry(self.root, textvariable = self.additional, width = 90, justify ='left')
        self.additional_text_box.grid(row = row, column = 1, padx = 5, sticky = 'EW', columnspan = 4)
        row += 1
        self.root.grid_rowconfigure(row,weight=1)
        Tkinter.Label(self.root, text = "Add Unconfirmed").grid(row = row, column = 0, padx = 15, pady = 5, sticky = 'E')
        self.add_unconfirmed_box.grid(row = row, column = 1, padx = 5, sticky = 'W')
        row += 1
        self.root.grid_rowconfigure(row,weight=1)
        Tkinter.Label(self.root, text = "Unconfirmed Planets").grid(row = row, column = 0, sticky = 'E', padx = 15, pady = 5)
        self.unconfirmed_text_box = Tkinter.Entry(self.root, textvariable = self.unconfirmed, width = 90, justify ='left')
        self.unconfirmed_text_box.grid(row = row, column = 1, padx = 5, sticky = 'EW', columnspan = 4)
        row += 1
        self.root.grid_rowconfigure(row,weight=1)
        Tkinter.Label(self.root, text = "Removed Planets").grid(row = row, column = 0, sticky = 'E', padx = 15, pady = 5)
        self.removed_text_box = Tkinter.Entry(self.root, textvariable = self.removed, width = 90, justify ='left')
        self.removed_text_box.grid(row = row, column = 1, padx = 5, sticky = 'EW',columnspan = 4)
        row += 1
        self.root.grid_rowconfigure(row,weight=1)
        self.run_button = Tkinter.Button(self.root, text =' Run Dynamite', command = self.run_dynamite)
        self.run_button.grid(row = row, column = 0, padx = 10, pady = 5, sticky = 'W')     
        self.plot_button = Tkinter.Button(self.root, text = 'Run Plots Only', command = self.plot_dynamite)
        self.plot_button.grid(row = row, column = 1, padx = 0, sticky = 'W')
        self.targets_editor_button = Tkinter.Button(self.root, text = 'Targets Editor', command = self.run_targets_editor)
        self.targets_editor_button.grid(row = row, column = 2, padx = 20, sticky = 'W')
        self.load_button = Tkinter.Button(self.root, text = 'Load New Config', command = self.load_config_file)
        self.load_button.grid(row = row, column = 3, padx = 0, sticky = 'W')
        self.edit_button = Tkinter.Button(self.root, text = 'Edit Config File', command = lambda: self.display_editor(self.root))
        self.edit_button.grid(row = row, column = 4, padx = 0, sticky = 'W')
        self.root.grid_rowconfigure(row, weight=1)
        self.exit_button = Tkinter.Button(self.root, text = 'Exit', command = self.exit_program)
        self.exit_button.grid(row = row, column = 5, padx = 5)
        
        self.root.grid_columnconfigure(0,weight=1)
        self.root.grid_columnconfigure(1,weight=1)
        self.root.grid_columnconfigure(2,weight=1)
        self.root.grid_columnconfigure(3,weight=1)
        self.root.grid_columnconfigure(4,weight=1)
        
        self.setup_widget_values()
        self.root.protocol("WM_DELETE_WINDOW", self.exit_program)
        self.root.mainloop()
    
    def set_plot_save_false(self, *args):
        self.save_plots_box.current(1)
        
    def check_valid_entries(self):
        """check that entries are valid lists"""
        for l in self.check_valid_entries_list:
            if not self.validate_entry(l.get()):
                return tkMB.askyesno(title="POSSIBLE INVALID ENTRY", message= "Possible Invalid Entry for " + self.find_config_key(l) + " -> " + l.get() + "\nDo you want to Continue?" )
        return True    
    
    def run_dynamite(self):
        """run dynamite program"""
        if self.config_file == None:
            tkMB.showerror("Missing Config File", "Cannot RUN Dynamite")
            return
        if self.run_entries_validator():
            if self.check_valid_entries():
                self.save_config_data()
                self.send_kill_signal()
                if len(self.run_editors) > 0:
                    if (tkMB.askyesno("Close Output Window(s)", "Do you want to close previous RUN output window(s)?")):
                        for e in self.run_editors:
                            e.destroy()
                        self.run_editors = {}
                self.editors = self.run_editors
                thread = threading.Thread(target=self.execute_program, args=("DYNAMITE", "dynamite_v2.py " + self.config_file))
                thread.daemon = True
                thread.start()
            return

    def run_entries_validator(self):
        """ Validate run entries"""
        if self.show_plots_box.get() == 'True':
            return self.plot_entries_validator()
        return True
    
    def plot_dynamite(self):
        """run plots program"""
        if self.config_file == None:
            tkMB.showerror("Missing Config File", "Cannot PLOT Dynamite")
            return
        if self.save_plots_box.get() == 'False':
            tkMB.showerror("CANNOT RUN PLOTS", "Use Saved Data is set to False")
            return
        if self.plot_entries_validator():
            if self.check_valid_entries():
                self.save_config_data()
                self.send_kill_signal()
                if len(self.plot_editors) > 0:
                    if (tkMB.askyesno("Close Output Window(s)", "Do you want to close previous PLOT output window(s)?")):
                        for e in self.plot_editors:
                            e.destroy()
                        self.plot_editors = {}
                self.editors = self.plot_editors
                thread = threading.Thread(target=self.execute_program, args=("DYNAMITE PLOTS", "dynamite_plots_v2.py " + self.config_file))
                thread.daemon = True
                thread.start()
        return

    def run_targets_editor(self):
        """run targets editor program"""
        dynamite_targets_editor.dynamite_targets_editor(self.root, "dynamite_targets_v2.py", None, None)

    def plot_entries_validator(self):
        """valid plot entries"""
        if self.mode_box.get() == "single":
            if self.plot_p_r_box.get() == "True":
                tkMB.showerror(title="INVALID ENTRY, PLOT ABORTED", message= " Mode must NOT be single when Plot P R is True")
                return False
            if self.plot_tdtp_box.get() == "True":
                tkMB.showerror(title="INVALID ENTRY, PLOT ABORTED", message= " Mode must NOT be single when Plot TD-TP is True")
                return False
            if self.plot_deltas_box.get() == "True" :
                tkMB.showerror(title="INVALID ENTRY, PLOT ABORTED", message= " Mode must NOT be single when Plot Deltas is True")
                return False
            if self.plot_ratios_box.get() == "True":
                tkMB.showerror(title="INVALID ENTRY, PLOT ABORTED", message= " Mode must NOT be single when Plot Ratios is True")
                return False
        else:
            if self.plot_hist_box.get() == "True":
                tkMB.showerror(title="INVALID ENTRY, PLOT ABORTED", message= " Mode must be single when Plot Hist is True")
                return False
            if self.plot_pdf_box.get() == "True":
                tkMB.showerror(title="INVALID ENTRY, PLOT ABORTED", message= " Mode must be single when Plot PDF is True")
                return False
            if self.show_rearth_box.get() == "True":
                tkMB.showerror(title="INVALID ENTRY, PLOT ABORTED", message= " Mode must be single when Show Rearth is True")
                return False
        return True
       
    def validate_key(self, value, text):
        """Validates that the keystrokes for the text boxes are allowed."""
        self.save_plots_box.current(1)
        if (text in '0123456789.'):
            try:
                if (len(value) > 0):
                    float(value)
                return True
            except ValueError:
                tkMB.showerror(title="Invalid Entry " + text, message="Entry \"" + str(value) + "\" must be a number")
                return False
        else:
            tkMB.showerror(title="Invalid Entry " + text, message="Entry \"" + str(value) + "\" must be a number")
            return False
        
    def validate_entry(self, value):
        """Validates the entry."""  
        try:   
            if isinstance(ast.literal_eval(value), list):
                return True    
        except:
            pass
        return False  
               
    def execute_program(self, title, command):
        """execute program"""
        editor = Tkinter.Toplevel(self.root)
        editor.bind('<ButtonRelease-3>', self.clicker, add='')
        text_pad = ScrolledText.ScrolledText(editor, width=100, height=40)
        editor.wm_title(title + "  -->  RUNNING")
        editor.protocol("WM_DELETE_WINDOW", lambda:self.close_editor_window(editor, title))
        text_pad.pack(fill=Tkinter.BOTH, expand=True)
        try:
            start_time = time.time()
            p = subprocess.Popen("python -u " + command, shell=True, cwd=".", stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            atexit.register(p.terminate)
            self.editors[editor] = p
            log_file = 'log_file.txt'
            if os.path.exists(log_file):
                log_file = 'log_file1.txt'
                if os.path.exists(log_file):
                    log_file = 'log_file.txt'
            with open(log_file, "w") as lfile:
                while self.root.winfo_exists() and editor in self.editors and self.editors[editor] != None:
                    out = p.stdout.readline()
                    if len(out) == 0:
                        break
                    try:
                        text_pad.insert(Tkconstants.END, out[0:2000])
                        if len(out) > 2000: 
                            text_pad.insert(Tkconstants.END, "\n\n!!!!OUTPUT TOO LARGE FOR DISPLAY. SEE FULL OUTPUT IN LOG_FILE.TXT!!!!\n\n")
                        text_pad.see(Tkconstants.END)
                        lfile.write(out.decode() + '\n')
                    except: pass
                    p.poll()
            self.editors[editor] = None
            editor.wm_title(title + "  -->  COMPLETED")
            message = "Program " + command + " has Completed in " + time.strftime('%H:%M:%S', time.gmtime(time.time() - start_time))
            text_pad.insert(Tkconstants.END, message)
            print(message)
        except:
            self.kill_processes(editor, title)
   
    def close_editor_window(self, editor, title):
        """Close process window and kill processes"""
        try:
            p = self.editors[editor]
            if p == None:
                editor.destroy()
                del self.editors[editor]
            else:
                if (tkMB.askyesno("RUNNING JOBS", "Do you want to terminate running jobs? This may take time!!!!")):
                    self.kill_processes(editor, title)
        except:
            editor.destroy()
            
    def kill_processes(self, editor, title):
        """kill processes"""
        self.send_kill_signal()
        p = self.editors[editor]
        if p == None:
            return
        try:
            p.terminate()
            p.kill()
        except Exception as e:
            print(e)
        self.editors[editor] = None
        editor.wm_title(title + "  -->  TERMINATED")
    
    def send_kill_signal(self):
        try:
            kill_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            kill_socket.connect(("127.0.0.1", 12345))
            kill_socket.close()
        except: pass
                          
    def load_config_data(self, filename):
        """Loads config file."""
        try:
            self.config_parameters = {}
            try:
                config_data = np.loadtxt(filename, dtype = str, delimiter = '::')
                with open(filename) as tfile:
                    data = tfile.readlines()
                self.tooltips = {}
                for l in data:
                    if l.startswith("#") and l.find(":") != -1 and not l.startswith("##"):
                        p = l.replace("'","")[2:-1].split(":")
                        self.tooltips[p[0]] = p[1]
            except IOError:
                tkMB.showerror("Missing Config File", "Cannot load config file error")
                return
            for i in range(len(config_data)):
                self.config_parameters[config_data[i, 0]] = config_data[i, 1]
            if len(self.config_parameters) == 0:
                raise TypeError("No Config entries loaded")
            self.config_file = filename
            self.root.title("DYNAMITE ---> " + self.config_file)
        except Exception as e:
            tkMB.showerror("Config File Load Failed", "Cannot load Config File " + filename)
            self.config_file = None
            self.root.title("DYNAMITE ---> NO CONFIG FILE LOADED")
            print(str(e))    
            print(self.config_parameters)
    
    def setup_widget_values(self):
        self.mc_chain.set(self.find_selected_item(self.find_config_key(self.mc_chain), self.mc_chain_values))
        self.radmin.set(self.find_selected_item(self.find_config_key(self.radmin), self.radmin_values))
        self.radmax.set(self.find_selected_item(self.find_config_key(self.radmax), self.radmax_values))
        self.additional.set(self.find_selected_item(self.find_config_key(self.additional), self.additional_values))
        self.unconfirmed.set(self.find_selected_item(self.find_config_key(self.unconfirmed), self.unconfirmed_values))
        self.removed.set(self.find_selected_item(self.find_config_key(self.removed), self.removed_values))
        self.plot_colors.set(self.find_selected_item(self.find_config_key(self.plot_colors), self.plot_colors_values))
        self.mode_box.current(self.find_selected_item(self.find_config_key(self.mode_box), self.mode_values))
        self.system_box.current(self.find_selected_item(self.find_config_key(self.system_box), self.system_values))
        self.period_box.current(self.find_selected_item(self.find_config_key(self.period_box), self.period_values))       
        self.radius_box.current(self.find_selected_item(self.find_config_key(self.radius_box), self.radius_values))       
        self.radtype_box.current(self.find_selected_item(self.find_config_key(self.radtype_box), self.radtype_values))   
        self.mass_radius_box.current(self.find_selected_item(self.find_config_key(self.mass_radius_box), self.mass_radius_values))   
        self.otegi_rho_box.current(self.find_selected_item(self.find_config_key(self.otegi_rho_box), self.otegi_rho_values))
        self.inclination_box.current(self.find_selected_item(self.find_config_key(self.inclination_box), self.inclination_values))    
        self.eccentricity_box.current(self.find_selected_item(self.find_config_key(self.eccentricity_box), self.eccentricity_values))    
        self.stability_box.current(self.find_selected_item(self.find_config_key(self.stability_box), self.stability_values))   
        self.plot_box.current(self.find_selected_item(self.find_config_key(self.plot_box), self.plot_values))
        self.show_plots_box.current(self.find_selected_item(self.find_config_key(self.show_plots_box), self.show_plots_values))
        self.save_plots_box.current(self.find_selected_item(self.find_config_key(self.save_plots_box), self.save_plots_values))
        self.plot_p_r_box.current(self.find_selected_item(self.find_config_key(self.plot_p_r_box), self.plot_p_r_values))
        self.plot_p_scale_box.current(self.find_selected_item(self.find_config_key(self.plot_p_scale_box), self.plot_p_scale_values))
        self.plot_tdtp_box.current(self.find_selected_item(self.find_config_key(self.plot_tdtp_box), self.plot_tdtp_values))
        self.plot_deltas_box.current(self.find_selected_item(self.find_config_key(self.plot_deltas_box), self.plot_deltas_values))
        self.plot_ratios_box.current(self.find_selected_item(self.find_config_key(self.plot_ratios_box), self.plot_ratios_values))
        self.plot_hist_box.current(self.find_selected_item(self.find_config_key(self.plot_hist_box), self.plot_hist_values))
        self.plot_pdf_box.current(self.find_selected_item(self.find_config_key(self.plot_pdf_box), self.plot_pdf_values))
        self.show_rearth_box.current(self.find_selected_item(self.find_config_key(self.show_rearth_box), self.show_rearth_values))
        self.use_mass_box.current(self.find_selected_item(self.find_config_key(self.use_mass_box), self.use_mass_values))
        self.ind_p_box.current(self.find_selected_item(self.find_config_key(self.ind_p_box), self.ind_p_values))
        self.ind_r_box.current(self.find_selected_item(self.find_config_key(self.ind_r_box), self.ind_r_values))
        self.ind_i_box.current(self.find_selected_item(self.find_config_key(self.ind_i_box), self.ind_i_values))
        self.add_unconfirmed_box.current(self.find_selected_item(self.find_config_key(self.add_unconfirmed_box), self.add_unconfirmed_values))
        CreateToolTip(self.mc_chain_text_box, self.root, self.tooltips, self.config_entries, self.mc_chain)
        CreateToolTip(self.mode_box, self.root, self.tooltips, self.config_entries)
        CreateToolTip(self.system_box, self.root, self.tooltips, self.config_entries)
        CreateToolTip(self.period_box, self.root, self.tooltips, self.config_entries)
        CreateToolTip(self.radius_box, self.root, self.tooltips, self.config_entries)
        CreateToolTip(self.radtype_box, self.root, self.tooltips, self.config_entries) 
        CreateToolTip(self.radmin_text_box, self.root, self.tooltips, self.config_entries, self.radmin) 
        CreateToolTip(self.radmax_text_box, self.root, self.tooltips, self.config_entries, self.radmax)
        CreateToolTip(self.mass_radius_box, self.root, self.tooltips, self.config_entries) 
        CreateToolTip(self.otegi_rho_box, self.root, self.tooltips, self.config_entries)
        CreateToolTip(self.inclination_box, self.root, self.tooltips, self.config_entries)
        CreateToolTip(self.eccentricity_box, self.root, self.tooltips, self.config_entries)
        CreateToolTip(self.stability_box, self.root, self.tooltips, self.config_entries)
        CreateToolTip(self.plot_box, self.root, self.tooltips, self.config_entries)
        CreateToolTip(self.show_plots_box, self.root, self.tooltips, self.config_entries)
        CreateToolTip(self.save_plots_box, self.root, self.tooltips, self.config_entries)
        CreateToolTip(self.plot_p_r_box, self.root, self.tooltips, self.config_entries)
        CreateToolTip(self.plot_p_scale_box, self.root, self.tooltips, self.config_entries)
        CreateToolTip(self.plot_tdtp_box, self.root, self.tooltips, self.config_entries)
        CreateToolTip(self.plot_deltas_box, self.root, self.tooltips, self.config_entries)
        CreateToolTip(self.plot_ratios_box, self.root, self.tooltips, self.config_entries)
        CreateToolTip(self.plot_hist_box, self.root, self.tooltips, self.config_entries)
        CreateToolTip(self.plot_pdf_box, self.root, self.tooltips, self.config_entries)
        CreateToolTip(self.show_rearth_box, self.root, self.tooltips, self.config_entries)
        CreateToolTip(self.use_mass_box, self.root, self.tooltips, self.config_entries)
        CreateToolTip(self.ind_p_box, self.root, self.tooltips, self.config_entries)
        CreateToolTip(self.ind_r_box, self.root, self.tooltips, self.config_entries)
        CreateToolTip(self.ind_i_box, self.root, self.tooltips, self.config_entries)
        CreateToolTip(self.add_unconfirmed_box, self.root, self.tooltips, self.config_entries)
        CreateToolTip(self.plot_colors_text_box, self.root, self.tooltips, self.config_entries, self.plot_colors)
        CreateToolTip(self.additional_text_box, self.root, self.tooltips, self.config_entries, self.additional)
        CreateToolTip(self.unconfirmed_text_box, self.root, self.tooltips, self.config_entries, self.unconfirmed)
        CreateToolTip(self.removed_text_box, self.root, self.tooltips, self.config_entries, self.removed)
        CreateToolTip(self.run_button, self.root, "Executes the dynamite.py script", None)
        CreateToolTip(self.plot_button, self.root, "Executes the dynamite_plots.py script", None)
        CreateToolTip(self.load_button, self.root, "Loads a new dynamite config file into the GUI", None)
        CreateToolTip(self.edit_button, self.root, "Displays an editor for the currently loaded dynamite config file", None)
        CreateToolTip(self.targets_editor_button, self.root, "Displays an editor for the targets dictionary in dynamite_targets.py", None)
        CreateToolTip(self.exit_button, self.root, "Exits the GUI and stops all python sub processes", None)
        
    def find_selected_item(self, item,  values):
        """Find selected combo box entry"""
        try:
            item = self.config_parameters[item]
        except:
            return 0 if isinstance(values, list) else values
        if not isinstance(values, list):
            return item
        for i in range(len(values)):
            if str(values[i]) == str(item):
                return i
        return 0
    
    def find_config_key(self, widget):
        """Find Config Key in Config Parameters"""
        for (k, v) in self.config_entries.items():
            if v == widget:
                return str(k)
        return None
        
    def save_config_data(self):
        """Saves config file."""
        try:
            data = self.diff_config_file()
            if data == None:
                print("No updates to config file")
                return
            print("Updates to config file")   
            with open(self.config_file, "w") as cfile:
                file_config_entries = []
                for l in data:
                    if "::" not in l:
                        cfile.write(l)
                    else:
                        parts = l.strip().split("::")
                        file_config_entries.append(parts[0])
                        found = False
                        for entry in self.config_entries.keys():
                            if parts[0] == entry:
                                found = True
                                cfile.write(parts[0] + "::" + str(self.config_entries[entry].get()) + "\n")
                                break
                        if not found:
                            cfile.write(l)
                for entry in self.config_entries.keys():
                    found = False
                    for c_entry in file_config_entries:
                        if entry == c_entry:
                            found = True
                            break
                    if not found:
                        cfile.write(entry + "::" + str(self.config_entries[entry].get()) + "\n")
        except Exception as e:
            print(str(e))    

    def diff_config_file(self):
        """ Check if there are config changes"""
        try:
            with open(self.config_file, "r") as cfile:
                data = cfile.readlines()
            changed = False
            c_count = 0
            for l in data:
                if "::" in l:
                    parts = l.strip().split("::")
                    c_count += 1
                    for entry in self.config_entries.keys(): 
                        if entry == "plt_colors":
                            self.config_entries[entry].set(self.config_entries[entry].get().replace("#",''))
                        if parts[0] == entry and str(parts[1]) != str(self.config_entries[entry].get()):
                            changed = True
            if c_count != len(self.config_entries):
                tkMB.showwarning(title="Number of Config entries differ", message= "GUI config Items " + str(len(self.config_entries)) + " Config File Entries " + str(c_count))
            return data if changed or c_count != len(self.config_entries) else None
        except Exception as e:
            print(str(e))
            return None
                       
    def exit_program(self):
        """Exits program.""" 
        self.send_kill_signal()
        sys.exit()
    
    def load_config_file(self):
        """Loads config file from disk."""
        file_handler = tkFD.askopenfilename()       
        if (file_handler == ""):
            return
        self.load_config_data(file_handler)
        self.setup_widget_values()
        
    def save_command(self, editor, text_pad):
        """Saves the file in the editor."""
        options = {}
        options['initialfile'] = self.config_file
        file_handler = tkFD.asksaveasfile(mode = 'w', **options)
        
        if (file_handler != None):
            self.config_file = file_handler.name
            file_handler.write(text_pad.get('1.0', Tkinter.END + '-1c').replace('\r\n', '\n'))
            file_handler.close()
            editor.destroy()
            self.load_config_data(file_handler)
            self.setup_widget_values()
            
    def exit_command(self, editor, text_pad):
        """Exits the config editor."""
        if (text_pad.get('1.0', Tkinter.END + '-1c') != self.contents):
            if (tkMB.askyesno("Save Changes", "Do you want to save changes before exit?")):
                self.save_command(editor, text_pad)
        editor.destroy()
           
    def display_editor(self, root):
        """Displays editor."""
        editor = Tkinter.Toplevel(root)
        editor.bind('<ButtonRelease-3>', self.clicker, add = '')
        text_pad = ScrolledText.ScrolledText(editor, width = 120, height = 45)
        editor.wm_title("Config File Editor for " + self.config_file)
        editor.protocol("WM_DELETE_WINDOW", lambda: self.exit_command(editor, text_pad))
        menu = Tkinter.Menu(editor)
        editor.config(menu = menu)
        filemenu = Tkinter.Menu(menu, tearoff = 0)
        filemenu.add_command(label = "Save", command = lambda: self.save_command(editor, text_pad))
        filemenu.add_command(label = "Exit", command = lambda: self.exit_command(editor, text_pad))
        menu.add_cascade(label = "File", menu = filemenu)
        editmenu = Tkinter.Menu(None, tearoff =0, takefocus = False) 
        editmenu.add_command(label = 'Cut', command = (lambda: self.click_cut(editor))) 
        editmenu.add_command(label = 'Copy', command = (lambda: self.click_copy(editor))) 
        editmenu.add_command(label = 'Paste', command = (lambda: self.click_paste(editor))) 
        editmenu.add_command(label = 'Select All', command = (lambda: self.click_select_all(editor))) 
        menu.add_cascade(label = "Edit", menu = editmenu)
        text_pad.pack(fill=ScrolledText.BOTH, side= ScrolledText.LEFT, expand=True)
        file_handler = open(self.config_file)
        self.contents = file_handler.read().replace('\r\n', '\n')
        text_pad.insert('1.0', self.contents)
        file_handler.close()
        text_pad.focus_set()

    def click_copy(self, e, apnd = False): 
        e.event_generate('<Control-c>')

    def click_cut(self, e):
        e.event_generate('<Control-x>') 

    def click_paste(self,e):
        e.event_generate('<Control-v>')

    def click_select_all(self,e):
        e.event_generate('<Control-a>')
       
    def clicker(self, e):
        """Shows events on right-mouse clicks for editor."""
        rmenu = Tkinter.Menu(None, tearoff = 0, takefocus = False) 
        rmenu.add_command(label = 'Cut', command = (lambda: self.click_cut(e.widget))) 
        rmenu.add_command(label = 'Copy', command = (lambda: self.click_copy(e.widget))) 
        rmenu.add_command(label = 'Paste', command = (lambda: self.click_paste(e.widget))) 
        rmenu.add_command(label = 'Select All', command = (lambda: self.click_select_all(e.widget))) 
        rmenu.tk_popup(e.x_root + 40, e.y_root + 10, entry = "0")

if __name__ == '__main__':      
    DynamiteGUI()

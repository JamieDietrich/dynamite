### DYNAmical Multi-planet Injection TEster (DYNAMITE) ###
### GUI ###
### Jeremy Dietrich ###
### jdietrich1@email.arizona.edu ###
### 2020 July 20 ###
### Version 1.2 ###
### Dietrich & Apai (2020), Astronomical Journal in press ###
### http://arxiv.org/pdf/2007.06745.pdf ###

import os
import sys
import ast
import shlex
import atexit
import threading
import subprocess
import scipy as sp
import tkinter as Tkinter
from tkinter import ttk
from tkinter import filedialog as tkFD
from tkinter import messagebox as tkMB
import tkinter.scrolledtext as ScrolledText
from tkinter import constants as Tkconstants

class DynamiteGUI:
    
    def __init__(self):
        """Creates an instance of the GUI with instance variables."""
        try:
            from dynamite_targets import dynamite_targets
            self.targets = dynamite_targets()
        except:
            self.targets = None
            
        self.config_entries = {}
        self.editor_entries = {}

        self.mode_values = ["single", "TESS", "Kepler", "all", "test"]
        self.period_values = ["epos", "syssim"]
        self.radius_values = ["epos", "syssim"]
        self.radtype_values = ["powerlaw", "clustered"]
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
        self.ind_p_values = ["linear_zoom", "linear", "log"]
        self.ind_r_values = ["linear_zoom", "linear"]
        self.ind_i_values = ["full", "truncated"]
                       
        self.root = Tkinter.Tk()
        
        ttk.Style().map("TCombobox", fieldbackground = [('!readonly', '#FFFFFF'),('readonly', '#FFFFFF'),],)
        
        if (len(sys.argv)  >=  2):
            self.config_file = (sys.argv[1])
        else:
            self.config_file = "dynamite_config.txt"
        self.root.title("DYNAMITE ---> " + self.config_file)
        
        self.mc_chain = Tkinter.StringVar()
        self.mc_chain.set(10000)
        self.config_entries["MC_chain"] = self.mc_chain
        self.system = Tkinter.StringVar()
        self.system.set("")
        self.config_entries["system"] = self.system
        self.radmin = Tkinter.StringVar()
        self.radmin.set(0.1)
        self.config_entries["radmin"] = self.radmin
        self.radmax = Tkinter.StringVar()
        self.radmax.set(5)
        self.config_entries["radmax"] = self.radmax
        self.additional = Tkinter.StringVar()
        self.additional.set("[[]]")
        self.config_entries["additional"] = self.additional
        self.unconfirmed = Tkinter.StringVar()
        self.unconfirmed.set("[[]]")
        self.config_entries["unconfirmed"] = self.unconfirmed
        self.removed = Tkinter.StringVar()
        self.removed.set("[]")
        self.config_entries["removed"] = self.removed
        self.plot_colors = Tkinter.StringVar()
        self.plot_colors.set('["377eb8", "d95f02", "4daf4a", "984ea3", "e41a1c", "e6ab02"]')
        self.config_entries["plt_colors"] = self.plot_colors

        self.mode_box = ttk.Combobox(self.root, width = 10, justify = Tkconstants.RIGHT)
        self.mode_box["values"] = self.mode_values
        self.mode_box.current(0)
        self.config_entries["mode"] = self.mode_box
        self.mode_box.bind('<<ComboboxSelected>>')
        self.mode_box.configure(state = 'readonly')
        
        self.period_box = ttk.Combobox(self.root, width = 10, justify = Tkconstants.RIGHT)
        self.period_box["values"] = self.period_values
        self.period_box.current(0)
        self.config_entries["period"] = self.period_box
        self.period_box.bind('<<ComboboxSelected>>')
        self.period_box.configure(state = 'readonly')
        
        self.radius_box = ttk.Combobox(self.root, width = 10, justify = Tkconstants.RIGHT)
        self.radius_box["values"] = self.radius_values
        self.radius_box.current(1)
        self.config_entries["radius"] = self.radius_box
        self.radius_box.bind('<<ComboboxSelected>>')
        self.radius_box.configure(state = 'readonly')
        
        self.radtype_box = ttk.Combobox(self.root, width = 10, justify = Tkconstants.RIGHT)
        self.radtype_box["values"] = self.radtype_values
        self.radtype_box.current(1)
        self.config_entries["radtype"] = self.radtype_box
        self.radtype_box.bind('<<ComboboxSelected>>')
        self.radtype_box.configure(state = 'readonly')

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
        
        self.ind_p_box = ttk.Combobox(self.root, width = 10, justify = Tkconstants.RIGHT)
        self.ind_p_box["values"] = self.ind_p_values
        self.ind_p_box.current(2)
        self.config_entries["ind_P"] = self.ind_p_box
        self.ind_p_box.bind('<<ComboboxSelected>>')
        self.ind_p_box.configure(state = 'readonly')
        
        self.ind_r_box = ttk.Combobox(self.root, width = 10, justify = Tkconstants.RIGHT)
        self.ind_r_box["values"] = self.ind_r_values
        self.ind_r_box.current(1)
        self.config_entries["ind_R"] = self.ind_r_box
        self.ind_r_box.bind('<<ComboboxSelected>>')
        self.ind_r_box.configure(state = 'readonly')
        
        self.ind_i_box = ttk.Combobox(self.root, width = 10, justify = Tkconstants.RIGHT)
        self.ind_i_box["values"] = self.ind_i_values
        self.ind_i_box.current(0)
        self.config_entries["ind_i"] = self.ind_i_box
        self.ind_i_box.bind('<<ComboboxSelected>>')
        self.ind_i_box.configure(state = 'readonly')
        
        vcmd = (self.root.register(self.validate_key), '%P', '%S')
        vcmd1 = (self.root.register(self.validate_system), '%P', '%S')
        vcmd2 = (self.root.register(self.remove_hash_symbol), '%P', '%S')
        
        self.check_valid_entries_list = [self.additional, self.unconfirmed, self.removed, self.plot_colors]
        
        row = 0
        Tkinter.Label(self.root, text = "V 1", font = ("Helvetica", 7)).grid(row = row, column = 4, sticky = 'EN')
        Tkinter.Label(self.root, text = "PREDICTION PARAMETERS").grid(row = row, column = 1, sticky = 'E', padx = 60, pady = 5, columnspan = 2)
        row += 1
        Tkinter.Label(self.root, text = "MC Chain").grid(row = row, column = 0, padx = 5, sticky = 'W')
        self.mc_chain_text_box = Tkinter.Entry(self.root, textvariable = self.mc_chain, width = 10, validate = "key", validatecommand = vcmd, justify ='center')
        self.mc_chain_text_box.grid(row = row, column = 0, padx = 0, sticky = 'E')
        Tkinter.Label(self.root, text = "Mode").grid(row = row, column = 1, padx = 20, sticky = 'W')
        self.mode_box.grid(row = row, column = 1, padx = 0, sticky = 'E')
        Tkinter.Label(self.root, text = "System").grid(row = row, column = 2, padx = 10, sticky ='W')
        self.system_text_box = Tkinter.Entry(self.root, textvariable = self.system, width = 10, validate = "focusout", validatecommand = vcmd1, justify ='center')
        self.system_text_box.grid(row = row, column = 2, padx = 0, sticky = 'E')
        Tkinter.Label(self.root, text = "Period").grid(row = row, column = 3, padx = 10, pady = 5, sticky = 'W')
        self.period_box.grid(row = row, column = 3, padx = 20, sticky = 'E')
        row += 1
        Tkinter.Label(self.root, text = "Radius").grid(row = row, column = 0, padx = 5, sticky = 'W', pady = 5)
        self.radius_box.grid(row = row, column = 0, padx = 0, sticky = 'E')
        Tkinter.Label(self.root, text = "Rad Type").grid(row = row, column = 1, padx = 20, sticky ='W')
        self.radtype_box.grid(row = row, column = 1, sticky = 'E')
        Tkinter.Label(self.root, text = "Rad Min").grid(row = row, column = 2, sticky = 'W', padx = 10)
        self.radmin_text_box = Tkinter.Entry(self.root, textvariable = self.radmin, width = 8, validate = "key", validatecommand = vcmd, justify ='center')
        self.radmin_text_box.grid(row = row, column = 2, padx = 20, sticky = 'E')
        Tkinter.Label(self.root, text = "Rad Max").grid(row = row, column = 3, sticky = 'W', padx = 10)
        self.radmax_text_box = Tkinter.Entry(self.root, textvariable = self.radmax, width = 8, validate = "key", validatecommand = vcmd, justify ='center')        
        self.radmax_text_box.grid(row = row, column = 3, padx = 20, sticky = 'E')
        row += 1
        Tkinter.Label(self.root, text = "PLOTTING PARAMETERS").grid(row = row, column = 1, sticky = 'E', padx = 60, pady = 5, columnspan = 2)
        row += 1
        Tkinter.Label(self.root, text = "Create Plots").grid(row = row, column = 0, padx = 5, sticky = 'W', pady = 5)
        self.plot_box.grid(row = row, column = 0, padx = 10, sticky = 'E')
        Tkinter.Label(self.root, text = "Show Plots").grid(row = row, column = 1, padx = 0, sticky = 'W')
        self.show_plots_box.grid(row = row, column = 1, padx = 0, sticky = 'E')
        Tkinter.Label(self.root, text = "Plot P R").grid(row = row, column = 2, sticky = 'W', padx = 10)
        self.plot_p_r_box.grid(row = row, column = 2, padx = 0, sticky = 'E')
        Tkinter.Label(self.root, text = "Use Saved Data").grid(row = row, column = 3, padx = 5, sticky = 'W')
        self.save_plots_box.grid(row = row, column = 3, padx = 0, sticky = 'E')
        row += 1
        Tkinter.Label(self.root, text = "Plot P Scale").grid(row = row, column = 0, padx = 5, sticky = 'W', pady = 5)
        self.plot_p_scale_box.grid(row = row, column = 0, padx = 10, sticky = 'E')
        Tkinter.Label(self.root, text = "Plot P(t)").grid(row = row, column = 1, padx = 10, sticky = 'W')
        self.plot_tdtp_box.grid(row = row, column = 1, padx = 0, sticky = 'E')
        Tkinter.Label(self.root, text = "Plot Deltas").grid(row = row, column = 2, padx = 10, sticky = 'W')
        self.plot_deltas_box.grid(row = row, column = 2, padx =  0, sticky = 'E')
        Tkinter.Label(self.root, text = "Plot Ratios").grid(row = row, column = 3, sticky = 'W', padx = 10)
        self.plot_ratios_box.grid(row = row, column = 3, padx = 0, sticky = 'E')
        row += 1
        Tkinter.Label(self.root, text = "Plot Hist").grid(row = row, column = 0, padx = 5, sticky = 'W', pady =5)
        self.plot_hist_box.grid(row = row, column = 0, padx = 10, sticky = 'E')
        Tkinter.Label(self.root, text = "Plot PDF").grid(row = row, column = 1, padx = 10, sticky = 'W')
        self.plot_pdf_box.grid(row = row, column = 1, padx = 0, sticky = 'E')
        Tkinter.Label(self.root, text = "Show Rearth").grid(row = row, column = 2, padx = 10, sticky = 'W')
        self.show_rearth_box.grid(row = row, column = 2, padx =0, sticky = 'E')
        row += 1
        Tkinter.Label(self.root, text = "IND P").grid(row = row, column = 0, sticky = 'W', padx = 5, pady = 5)
        self.ind_p_box.grid(row = row, column = 0, padx = 10, sticky = 'E')
        Tkinter.Label(self.root, text = "IND R").grid(row = row, column = 1, padx = 10, sticky = 'W')
        self.ind_r_box.grid(row = row, column = 1, padx = 0, sticky = 'E')
        Tkinter.Label(self.root, text = "IND I").grid(row = row, column = 2, padx = 10, sticky = 'W')
        self.ind_i_box.grid(row = row, column = 2, padx = 0, sticky = 'E')
        row += 1
        Tkinter.Label(self.root, text = "Plot Colors").grid(row = row, column = 0, sticky = 'W', padx = 5, pady = 5)
        self.plot_colors_text_box = Tkinter.Entry(self.root, textvariable = self.plot_colors, width = 70, validate = "focusout", validatecommand = vcmd2, justify ='left')
        self.plot_colors_text_box.grid(row = row, column = 0, padx = 40, sticky = 'E',columnspan = 4)
        row += 1
        Tkinter.Label(self.root, text = "SYSTEM ARCHITECTURE PARAMETERS").grid(row = row, column = 1, sticky = 'E', padx = 0, pady = 5, columnspan = 2)
        row += 1
        Tkinter.Label(self.root, text = "Additional Planets").grid(row = row, column = 0, sticky = 'W', padx = 5, pady = 5)
        self.additional_text_box = Tkinter.Entry(self.root, textvariable = self.additional, width = 74, justify ='left')
        self.additional_text_box.grid(row = row, column = 0, padx = 5, sticky = 'E', columnspan = 5)
        row += 1
        Tkinter.Label(self.root, text = "Unconfirmed Planets").grid(row = row, column = 0, sticky = 'W', padx = 5, pady = 5)
        self.unconfirmed_text_box = Tkinter.Entry(self.root, textvariable = self.unconfirmed, width = 74, justify ='left')
        self.unconfirmed_text_box.grid(row = row, column = 0, padx = 5, sticky = 'E', columnspan = 5)
        row += 1
        Tkinter.Label(self.root, text = "Removed Planets").grid(row = row, column = 0, sticky = 'W', padx = 5, pady = 5)
        self.removed_text_box = Tkinter.Entry(self.root, textvariable = self.removed, width = 74, justify ='left')
        self.removed_text_box.grid(row = row, column = 0, padx = 5, sticky = 'E',columnspan = 5)
        row += 1
        Tkinter.Button(self.root, text =' Run Dynamite', command = self.run_dynamite).grid(row = row, column = 0, padx = 30, pady = 5)
        Tkinter.Button(self.root, text = 'Run Plots Only', command = self.run_plots).grid(row = row, column = 1, padx = 30)
        Tkinter.Button(self.root, text = 'Load New Config', command = self.load_config_file).grid(row = row, column = 2, padx = 30)
        Tkinter.Button(self.root, text = 'Edit Config File', command = lambda: self.display_editor(self.root)).grid(row = row, column = 3, padx = 30)
        Tkinter.Button(self.root, text = 'Exit', command = self.exit_program).grid(row = row, column = 4, padx = 10)
       
        self.load_config_data()         
        self.root.resizable(0, 0)
        self.root.mainloop()
           
    def check_valid_entries(self):
        """check that entries are valid lists"""
        for l in self.check_valid_entries_list:
            if not self.validate_entry(l.get()):
                tkMB.showerror(title="Invalid Entry for " + self.find_config_key(l) + " cannot Run Dynamite", message=l.get() )
                return False
        return True    

    def run_dynamite(self):
        """run dynamitr program"""
        if self.check_valid_entries():
            self.save_config_data()
            thread = threading.Thread(target=self.execute_program, args=("DYNAMITE", "dynamite.py"))
            thread.daemon = True
            thread.start()
        return

    def run_plots(self):
        """run plots program"""
        if self.check_valid_entries():
            thread = threading.Thread(target=self.execute_program, args=("DYNAMITE PLOTS", "dynamite_plots.py"))
            thread.daemon = True
            thread.start()
            self.save_config_data()
        return
       
    def validate_key(self, value, text):
        """Validates that the keystrokes for the text boxes are allowed."""        
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
            return isinstance(ast.literal_eval(value), list)
        except:
            return False  
        
    def validate_system(self, value, text):
        """Validates the system."""        
        if len(value) == 0:
            return True
        try:
            if self.targets == None or len(self.targets.get_targets("single", value, None, [])) > 0:
                return True
        except:
            pass
        tkMB.showerror(title="Invalid System Target", message="System " + value + " does NOT exist in targets dictionary")
        return False
    
    def remove_hash_symbol(self, value, text):
        """Remove # symbol from color string""" 
        self.plot_colors.set(value.replace("#", ""))
        return True
    
    def execute_program(self, title, program):
        """execute program"""
        editor = Tkinter.Toplevel(self.root)
        editor.bind('<ButtonRelease-3>', self.clicker, add='')
        text_pad = ScrolledText.ScrolledText(editor, width=140, height=40)
        editor.wm_title(title)
        editor.protocol("WM_DELETE_WINDOW", lambda:self.close_editor_window(editor))
        text_pad.pack(fill=Tkinter.BOTH, expand=True)
        try:
            p = subprocess.Popen(shlex.split("python -u " + program + " " + self.config_file, posix=False), shell=False, cwd=".", stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            atexit.register(p.terminate)
            self.editor_entries[editor] = p
            while self.root.winfo_exists() and editor and p != None:
                out = p.stdout.readline()
                if out == '' and p.poll() != None:
                    break
                text_pad.insert(Tkconstants.END, out)
                text_pad.see(Tkconstants.END)
        except Exception as e:
            self.kill_processes(p, editor)
            print(e)
   
    def close_editor_window(self, editor):
        """Close process window and kill processes"""
        p = self.editor_entries[editor]
        if p == None:
            editor.destroy()
            del self.editor_entries[editor]
        else:
            self.kill_processes(p, editor)
            if (tkMB.askyesno("Close Output Window", "Do you want to close output window?")):
                editor.destroy()
                del self.editor_entries[editor]
    
    def kill_processes(self, p, editor):
        """kill processes"""
        if p == None:
            return
        try:
            p.terminate()
            p.kill()
            os.kill(p.pid, 0)
        except Exception:
            pass
        self.editor_entries[editor] = None
                  
    def load_config_data(self):
        """Loads config file."""
        try:
            self.config_parameters = {}
            try:
                config_data = sp.loadtxt(self.config_file, dtype = str, delimiter = '::')
            except IOError:
                tkMB.showerror("Missing Config File", "Cannot load config file fatal error")
                sys.exit()
            for i in range(len(config_data)):
                self.config_parameters[config_data[i, 0]] = config_data[i, 1]

            self.mc_chain.set(self.config_parameters[self.find_config_key(self.mc_chain)])
            self.system.set(self.config_parameters[self.find_config_key(self.system)])
            self.radmin.set(self.config_parameters[self.find_config_key(self.radmin)])
            self.radmax.set(self.config_parameters[self.find_config_key(self.radmax)])
            self.additional.set(self.config_parameters[self.find_config_key(self.additional)])
            self.unconfirmed.set(self.config_parameters[self.find_config_key(self.unconfirmed)])
            self.removed.set(self.config_parameters[self.find_config_key(self.removed)])
            self.plot_colors.set(self.config_parameters[self.find_config_key(self.plot_colors)])
            self.mode_box.current(self.find_selected_item(self.config_parameters[self.find_config_key(self.mode_box)], self.mode_values))
            self.period_box.current(self.find_selected_item(self.config_parameters[self.find_config_key(self.period_box)], self.period_values))       
            self.radius_box.current(self.find_selected_item(self.config_parameters[self.find_config_key(self.radius_box)], self.radius_values))       
            self.radtype_box.current(self.find_selected_item(self.config_parameters[self.find_config_key(self.radtype_box)], self.radtype_values))       
            self.plot_box.current(self.find_selected_item(self.config_parameters[self.find_config_key(self.plot_box)], self.plot_values))
            self.show_plots_box.current(self.find_selected_item(self.config_parameters[self.find_config_key(self.show_plots_box)], self.show_plots_values))
            self.save_plots_box.current(self.find_selected_item(self.config_parameters[self.find_config_key(self.save_plots_box)], self.save_plots_values))
            self.plot_p_r_box.current(self.find_selected_item(self.config_parameters[self.find_config_key(self.plot_p_r_box)], self.plot_p_r_values))
            self.plot_p_scale_box.current(self.find_selected_item(self.config_parameters[self.find_config_key(self.plot_p_scale_box)], self.plot_p_scale_values))
            self.plot_tdtp_box.current(self.find_selected_item(self.config_parameters[self.find_config_key(self.plot_tdtp_box)], self.plot_tdtp_values))
            self.plot_deltas_box.current(self.find_selected_item(self.config_parameters[self.find_config_key(self.plot_deltas_box)], self.plot_deltas_values))
            self.plot_ratios_box.current(self.find_selected_item(self.config_parameters[self.find_config_key(self.plot_ratios_box)], self.plot_ratios_values))
            self.plot_hist_box.current(self.find_selected_item(self.config_parameters[self.find_config_key(self.plot_hist_box)], self.plot_hist_values))
            self.plot_pdf_box.current(self.find_selected_item(self.config_parameters[self.find_config_key(self.plot_pdf_box)], self.plot_pdf_values))
            self.show_rearth_box.current(self.find_selected_item(self.config_parameters[self.find_config_key(self.show_rearth_box)], self.show_rearth_values))
            self.ind_p_box.current(self.find_selected_item(self.config_parameters[self.find_config_key(self.ind_p_box)], self.ind_p_values))
            self.ind_r_box.current(self.find_selected_item(self.config_parameters[self.find_config_key(self.ind_r_box)], self.ind_r_values))
            self.ind_i_box.current(self.find_selected_item(self.config_parameters[self.find_config_key(self.ind_i_box)], self.ind_i_values))
    
            self.root.title("DYNAMITE ---> " + self.config_file)
        except Exception as e:
            print(str(e))    
            print(self.config_parameters)

    def find_config_key(self, widget):
        """Find Config Key in Config Parameters"""
        for (k, v) in self.config_entries.items():
            if v == widget:
                return str(k)
            
    def save_config_data(self):
        """Saves config file."""
        try:
            data = self.diff_config_file()
            if data == None:
                print("No updates to config file")
                return
            print("Updates to config file")   
            with open(self.config_file, "w") as cfile:
                for l in data:
                    if "::" not in l:
                        cfile.write(l)
                    else:
                        parts = l.strip().split("::")
                        for entry in self.config_entries.keys():        
                            if parts[0] == entry:
                                cfile.write(parts[0] + "::" + str(self.config_entries[entry].get()) + "\n")
                                break
        except Exception as e:
            print(str(e))    

    def diff_config_file(self):
        """ Check if there are config changes"""
        try:
            with open(self.config_file, "r") as cfile:
                data = cfile.readlines()
            for l in data:
                if "::" in l:
                    parts = l.strip().split("::")
                    for entry in self.config_entries.keys():     
                        if parts[0] == entry and str(parts[1]) != str(self.config_entries[entry].get()):
                            return data
            return None
        except Exception as e:
            print(str(e))
            return None
        
    def find_selected_item(self, item, items):
        """Find selected combo box entry"""
        for i in range(len(items)):
            if str(items[i]) == str(item):
                return i
        return 0
                       
    def exit_program(self):
        """Exits program."""  
        sys.exit()
    
    def load_config_file(self):
        """Loads config file from disk."""
        file_handler = tkFD.askopenfilename()       
        if (file_handler == ""):
            return
        self.config_file = file_handler
        self.load_config_data()
        
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
            self.load_config_data()
         
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
        text_pad.pack(fill=ScrolledText.BOTH, side= ScrolledText.LEFT, expand=True)
        file_handler = open(self.config_file)
        self.contents = file_handler.read().replace('\r\n', '\n')
        text_pad.insert('1.0', self.contents)
        file_handler.close()
        
    def clicker(self, e):
        """Shows events on right-mouse clicks for editor."""
        def click_copy(e, apnd = False): 
            e.widget.event_generate('<Control-c>')

        def click_cut(e):
            e.widget.event_generate('<Control-x>') 

        def click_paste(e):
            e.widget.event_generate('<Control-v>')

        rmenu = Tkinter.Menu(None, tearoff = 0, takefocus = False) 
        rmenu.add_command(label = 'Cut', command = (lambda: click_cut(e))) 
        rmenu.add_command(label = 'Copy', command = (lambda: click_copy(e))) 
        rmenu.add_command(label = 'Paste', command = (lambda: click_paste(e))) 
        rmenu.tk_popup(e.x_root + 40, e.y_root + 10, entry = "0")
        
DynamiteGUI()

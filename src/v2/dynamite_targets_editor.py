### DYNAmical Multi-planet Injection TEster (DYNAMITE) ###
### Targets Editor ###
### Jeremy Dietrich ###
### jdietrich1@email.arizona.edu ###
### 2022 January 8 ###
### Version 2.0 ###
### Dietrich & Apai (2020), AJ, 160, 107D ###
### Dietrich & Apai (2021), AJ, 161, 17D ###
### Dietrich, Apai, & Malhotra (2022), accepted to AJ ###
### https://iopscience.iop.org/article/10.3847/1538-3881/aba61d ###

# The command line call for this function is "python [editor] [targets file]"

import re
import sys
import copy
import time
import threading
import traceback
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

class dynamite_targets_editor():
    def __init__(self, root, filename, header, text_variables):
        self.rootMain=root
        self.headers_info = header
        self.text_variables = text_variables
        self.headers_length = len(self.headers_info)
        self.planet_fields  = len([x for x in self.headers_info if x[2] == 'P'])
        self.stellar_fields = len([x for x in self.headers_info if x[2] == 'S'])
        self.planet_offset = self.stellar_fields + 1 
        self.headers = False
        self.max_rows = int(root.winfo_screenheight()/40)-4
        self.top_row = 0
        self.bottom_row = self.top_row + self.max_rows
        self.reverse = False
        self.find_values_list = []
        self.replace_values_list = []
        self.label_list = []
        self.undo_list = []
        self.close_program = False
        self.root = tk.Toplevel(self.rootMain)
        self.title = "Targets Dictionary"
        self.root.title(self.title)
        self.root.configure(background="Gray")
        self.root.columnconfigure(0, weight=1)
        self.vcmd = (self.root.register(self.modify_entry), '%W', '%P', '%S')
        self.vcmd1 = (self.root.register(self.enable_replace_buttons), '%P')
        self.vcmd2 = (self.root.register(self.enable_find_buttons), '%P')
        self.root.geometry('%dx%d+%d+%d' % (root.winfo_screenwidth(), (self.max_rows*44)+5, 0, 0))
        self.previous_widget = None
        self.draw_grid()
        self.loaded_filename = filename
        if self.loaded_filename != None:
            self.load_file(self.loaded_filename)
        self.entry_boxes[0][0][0].focus_set()
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)
        self.root.resizable(1,0)
        self.display_rows()
        self.root.focus_force()
        self.root.mainloop()

    def draw_grid(self):     
        frame1 = tk.Frame(self.root, bg="light grey", bd=2, relief=tk.GROOVE)
        frame1.grid(row=0, column=0, sticky="new")
        for i in range(self.headers_length):
            frame1.columnconfigure(i, weight=1)
            label = tk.Label(frame1, width = self.headers_info[i][4], text=self.headers_info[i][3], font=('Arial',10,))
            label.grid(row=0, column=i, pady=0, padx=2, sticky="WE")
            label.bind('<ButtonRelease-1>',  self.sort_column)
            self.label_list.append(label)
        frame2 = tk.Frame(self.root)
        frame2.grid(row=1, column=0, sticky=tk.NSEW)
        frame2.columnconfigure(0, weight=1)
        self.canvas = tk.Canvas(frame2, bg="Gray")
        self.canvas.grid(row=0, column=0, sticky =tk.NSEW)
        self.canvas.columnconfigure(0, weight=1)
        self.entry_frame = tk.Frame(self.canvas, bg="light grey", bd=0)
        self.entry_frame.grid(row=0, column=0, sticky =tk.NSEW)
        self.entry_boxes = []
        for row in range(self.max_rows + 20):
            entry_list = []
            for col in range(self.headers_length):
                entry = tk.Entry(self.entry_frame, width = self.headers_info[col][4], validate = "key", validatecommand = self.vcmd, justify = 'center', fg='black', bg = "SlateGray1", font=('Arial',10,), readonlybackground ='light grey', state='readonly')
                entry.grid(row=row, column=col, sticky='news', padx =2)
                entry.bind("<Prior>", lambda event: self.move_display("page up"))
                entry.bind("<Next>",  lambda event: self.move_display("page dn"))
                entry.bind("<Up>",    lambda event: self.move_display("scroll up"))
                entry.bind("<Down>",  lambda event: self.move_display("scroll dn"))
                entry.bind("<Home>",  lambda event: self.move_display("home"))
                entry.bind("<End>",   lambda event: self.move_display("end"))
                entry.bind("<Escape>",lambda event: self.remove_find_popup())
                entry.bind("<Control-f>", lambda event: self.create_find_popup())
                entry.bind("<F1>", lambda event: self.create_help_popup())
                entry.bind("<Control-z>", lambda event: self.do_undo())
                entry.bind("<Control-m>", lambda event: self.create_modify_header_popup())
                entry.bind("<Control-l>", lambda event: self.load_file(None, False))
                entry.bind("<Control-L>", lambda event: self.load_file(None, True))
                entry.bind("<Control-s>", lambda event: self.save_file())
                entry.bind("<Control-i>", lambda event: self.import_csv_file())
                entry.bind("<Control-Delete>", lambda event: self.delete_all_targets())
                entry.bind("<FocusIn>", self.focus)
                self.entry_frame.columnconfigure(col, weight=1)
                if col == 0 or col == self.planet_offset:
                    entry.bind('<ButtonRelease-3>', self.clicker_target if col == 0 else self.clicker_planet)
                entry_list.append((entry, row, col, 0, 0))
            self.entry_boxes.append(entry_list)

    def create_modify_header_popup(self):
        if not hasattr(self, "modify_header_top") or not self.modify_header_top.winfo_exists():
            self.header_index = 0
            self.type = tk.StringVar()
            self.name = tk.StringVar()
            self.size = tk.StringVar()
            self.def_data = tk.StringVar()
            self.def_data.set("?")
            self.header_label = tk.StringVar()
            self.total_header_fields = self.headers_length
            header = []
            for i in range(self.headers_length):
                for j in range(3):
                    header.append(str(self.headers_info[i][j+2]) + ("" if i == self.headers_length-1 and j == 2 else "|"))
            self.header_label.set(' '.join(map(str, header)))
            self.swap_start_index = -1
            self.display_header_section(0)
            self.modify_header_top = tk.Toplevel(self.root)
            self.modify_header_top.attributes('-topmost', 'true')
            self.modify_header_top.title('Modify Editor Header Fields                                       Use Up and Down Arrows to Scroll Through Header Fields')
            vcmd = (self.root.register(self.validate_type), '%P', '%S')
            vcmd1 = (self.root.register(self.validate_name), '%P', '%S')
            vcmd2 = (self.root.register(self.validate_size), '%P', '%S')
            tk.Label(self.modify_header_top, textvariable=self.header_label).grid(row=0, column=0, sticky='w', padx =10, pady = 10,columnspan=8)
            tk.Label(self.modify_header_top, text="Type").grid(row=1, column=0, sticky='e', padx =10, pady = 10)
            self.type_entry = tk.Entry(self.modify_header_top, textvariable = self.type, width = 5, justify ='center', validate = "key", validatecommand = vcmd)
            self.type_entry.grid(row=1, column=1, sticky='w', padx = 10, pady = 10)
            self.type_entry.bind("<Escape>", lambda event: self.remove_find_popup())
            self.type_entry.bind("<Up>", lambda event: self.display_header_section(self.header_index +1))
            self.type_entry.bind("<Down>",  lambda event: self.display_header_section(self.header_index -1))
            tk.Label(self.modify_header_top, text="Title").grid(row=1, column=2, sticky='e', padx = 10, pady = 10)
            self.name_entry = tk.Entry(self.modify_header_top, textvariable = self.name, width = 15, justify ='center', validate = "key", validatecommand = vcmd1)
            self.name_entry.grid(row=1, column=3, sticky='w', padx =10, pady = 10)
            self.name_entry.bind("<Escape>", lambda event: self.remove_find_popup())
            self.name_entry.bind("<Up>", lambda event: self.display_header_section(self.header_index +1))
            self.name_entry.bind("<Down>",  lambda event: self.display_header_section(self.header_index -1))
            tk.Label(self.modify_header_top, text="Size").grid(row=1, column=4, sticky='e', padx = 10, pady = 10)
            self.size_entry = tk.Entry(self.modify_header_top, textvariable = self.size, width = 5, justify ='center', validate = "key", validatecommand = vcmd2)
            self.size_entry.grid(row=1, column=5, sticky='w', padx =10, pady = 10)
            self.size_entry.bind("<Escape>", lambda event: self.remove_find_popup())
            self.size_entry.bind("<Up>", lambda event: self.display_header_section(self.header_index +1))
            self.size_entry.bind("<Down>",  lambda event: self.display_header_section(self.header_index -1))
            tk.Label(self.modify_header_top, text="Default Data").grid(row=1, column=6, sticky='e', padx = 10, pady = 10)
            self.data_entry = tk.Entry(self.modify_header_top, textvariable = self.def_data, width = 10, justify ='center')
            self.data_entry.grid(row=1, column=7, sticky='w', padx =10, pady = 10)
            self.data_entry.bind("<Escape>", lambda event: self.remove_find_popup())
            self.data_entry.bind("<Up>", lambda event: self.display_header_section(self.header_index +1))
            self.data_entry.bind("<Down>",  lambda event: self.display_header_section(self.header_index -1))
            self.update_button = tk.Button(self.modify_header_top, text='Update',  command = lambda: self.do_update_header())
            self.update_button.grid(row=3, column=0, sticky='w', padx = 10, pady = 0)
            self.swap_button = tk.Button(self.modify_header_top, text='Swap Header',  command = lambda: self.do_swap_header())
            self.swap_button.grid(row=3, column=2, sticky='w', padx = 10, pady = 0)
            self.insert_button = tk.Button(self.modify_header_top, text='Insert Header After',  command = lambda: self.do_insert_header())
            self.insert_button.grid(row=3, column=4, sticky='w', padx = 10, pady = 0)
            self.delete_button = tk.Button(self.modify_header_top, text='Delete Header',  command = lambda: self.do_delete_header())
            self.delete_button.grid(row=3, column=6, sticky='w', padx = 10, pady = 0)
            self.cancel_button = tk.Button(self.modify_header_top, text='Cancel', command=self.remove_modify_header_popup)
            self.cancel_button.grid(row=3, column=7, sticky='e', padx = 10, pady = 0)
            self.modify_header_top.geometry("+%d+%d" % (self.root.winfo_x() + 200, self.root.winfo_y() + 200))
            self.name_entry.focus_set()
       
    def validate_type(self, value, text):
        if len(text) > 0:
            if text in "TSP":
                header = list(self.header_label.get().split(" "))
                header[self.header_index * 3] = value + "|" 
                self.header_label.set(' '.join(map(str, header)).replace("'","").replace(")",""))
            else:
                return False
        return True
    
    def validate_name(self, value, text):
        if len(text) > 0:
            header = list(self.header_label.get().split(" ")) 
            header[(self.header_index * 3) +1] = value + "|"
            self.header_label.set(' '.join(map(str, header)).replace("'","").replace(")",""))
        return True
    
    def validate_size(self, value, text):
        if len(text) > 0:
            if text.isnumeric():
                header = list(self.header_label.get().split(" "))
                header[(self.header_index * 3) + 2] = value + "|"
                self.header_label.set(' '.join(map(str, header)).replace("'","").replace(")",""))
            else:
                return False
        return True
    
    def display_header_section(self, index):
        if index < 0:
            index = 0
        if index == self.headers_length:
            index = self.headers_length -1
        self.header_index = index
        header_info = list(self.header_label.get().split("| "))
        self.type.set(header_info[index * 3])
        self.name.set(header_info[index * 3 + 1])
        self.size.set(header_info[index * 3 + 2])
    
    def do_swap_header(self):
        if self.swap_start_index == -1:
            self.swap_start_index = self.header_index
        else:
            header_info = list(self.header_label.get().split("| "))
            if header_info[self.swap_start_index*3] != header_info[self.header_index*3]:
                messagebox.showerror(title="Swap Header Error", message= "Cannot swap headers they are not the same type")
                self.swap_start_index = -1
                return
            header_info = list(self.header_label.get().split("| "))
            h_type = header_info[self.swap_start_index*3]
            title = header_info[self.swap_start_index*3 +1]
            size = header_info[self.swap_start_index*3 +2]
            header_info[self.swap_start_index*3] = header_info[self.header_index*3]
            header_info[self.swap_start_index*3 +1] = header_info[self.header_index*3 +1]
            header_info[self.swap_start_index*3 +2] = header_info[self.header_index*3 +2]
            header_info[self.header_index*3] = h_type
            header_info[self.header_index*3 +1] = title
            header_info[self.header_index*3 +2] = size
            self.header_label.set('| '.join(map(str, header_info)).replace("'","").replace(")",""))
            self.display_header_section(self.header_index)
            if messagebox.askyesno("Swap Data", "Do you want to swap the column data also?"):
                if self.text_variables != None:
                    for i in range(len(self.text_variables)):
                        if self.header_index <= self.stellar_fields+1 and self.swap_start_index <= self.stellar_fields+1:
                            temp = self.text_variables[i][self.swap_start_index]
                            self.text_variables[i][self.swap_start_index] = self.text_variables[i][self.header_index]
                            self.text_variables[i][self.header_index] = temp
                        else:
                            inx = self.planet_offset + self.header_index - self.stellar_fields -1
                            inx1 = self.planet_offset + self.swap_start_index - self.stellar_fields-1
                            while inx < len(self.text_variables[i]):
                                temp = self.text_variables[i][inx1]
                                self.text_variables[i][inx1] = self.text_variables[i][inx]
                                self.text_variables[i][inx] = temp
                                inx += self.planet_fields
                                inx1 += self.planet_fields
                    self.swap_start_index = -1
                    self.do_update_header()
            self.swap_start_index = -1
                     
    def do_delete_header(self):
        self.swap_start_index = -1
        self.total_header_fields -= 1
        header_info = list(self.header_label.get().split("| "))
        header = []
        inx = 0
        while inx < len(header_info):
            if inx == self.header_index*3:
                inx += 3
            else:
                header.append(header_info[inx] + ("|" if inx < len(header_info)-1 else ""))
                inx +=1
        self.header_label.set(' '.join(map(str, header)).replace("'","").replace(")",""))
        self.display_header_section(self.header_index)
        if messagebox.askyesno("Delete Data", "Do you want to delete the column data also?"):
            if self.text_variables != None:
                for i in range(len(self.text_variables)):
                    if self.header_index <= self.stellar_fields+1:
                        self.text_variables[i].pop(self.header_index)
                    else:
                        inx = self.planet_offset + self.header_index - self.stellar_fields-1
                        while inx < len(self.text_variables[i]):
                            self.text_variables[i].pop(inx)
                            inx += self.planet_fields-1
                self.do_update_header()
                    
    def do_insert_header(self):
        self.swap_start_index = -1
        self.total_header_fields += 1
        header_info = list(self.header_label.get().split("| "))
        header = []
        self.header_index += 1
        for i in range(len(header_info)):
            if i == self.header_index*3:
                header.append(self.type.get()+ "|")
                header.append("?" + "|")
                header.append(self.size.get() + "|")
            header.append(header_info[i] + ("|" if i < len(header_info)-1 else ""))
        self.header_label.set(' '.join(map(str, header)).replace("'","").replace(")",""))
        self.display_header_section(self.header_index)
        if messagebox.askyesno("Insert Data", "Do you want to insert the default column data also?"):
            if self.text_variables != None:
                for i in range(len(self.text_variables)):
                    if self.header_index <= self.stellar_fields+1:
                        self.text_variables[i].insert(self.header_index, ('"' + self.def_data.get() + '"',))
                    else:
                        inx = self.planet_offset + self.header_index - self.stellar_fields
                        while inx < len(self.text_variables[i]):
                            self.text_variables[i].insert(inx-1, ('"' + self.def_data.get() + '"',))
                            inx += self.planet_fields+1
                self.do_update_header()
              
    def do_update_header(self):
        header_pre_post =[('', ': '),('([', ', '),('', ', '),('', '],'),(' [', ', '),('', ', '),('', '],'),('', '')]
        headers = list(self.header_label.get().split("| "))
        header_info = [] 
        inx = 0
        state = -1
        for _ in range(self.total_header_fields):
            if state == -1: state = 0
            else: 
                if state == 0  and headers[inx] == "S": state = 1
                else:
                    if state == 1 and headers[inx] == "S": state = 2
                    else:
                        if state == 2 and headers[inx] == "S" and headers[inx+3] != "P": state = 2
                        else:
                            if state == 2 and headers[inx+3] == "P": state = 3
                            else:
                                if state == 3 and headers[inx] == "P": state = 4
                                else:
                                    if state == 4 and headers[inx] == "P" and headers[inx+3] != "C": state = 5
                                    else:
                                        if state == 5 and headers[inx+3] == "C": state = 6
                                        else:
                                            if state == 6 and headers[inx] == "C": state = 7
            header_info.append((header_pre_post[state][0], header_pre_post[state][1], headers[inx], headers[inx+1], headers[inx+2]))
            inx += 3
        with open(sys.argv[0], "r") as pfile:
            plines = pfile.readlines()
        with open(sys.argv[0], "w") as pfile:
            for line in plines:
                if line.startswith("        headerxx="):
                    line = "        headerxx=" + str(header_info) +"\n"
                pfile.write(line)
        self.root.destroy()
        dynamite_targets_editor(self.rootMain, None, header_info, self.text_variables)
        
    def remove_modify_header_popup(self):
        try: self.modify_header_top.destroy()
        except: pass
        try: self.help_top.destroy()
        except: pass
        
    def delete_all_targets(self):
        if (messagebox.askyesno("Delete All Targets", "Do you want to delete all targets?")):
            self.text_variables = []
            self.display_rows()
        
    def import_csv_file(self):
        files = filedialog.askopenfiles(mode = 'r')
        if files == None or len(files) == 0: return
        self.validate_active = False
        total_imported = 0
        total_read = 0
        for file in files:
            if ".csv" in file.name:
                tr, ti = self.process_csv_file(file.readlines(), total_read, total_imported)
                total_read += tr
                total_imported += ti
            file.close()
        self.display_rows()
        self.validate_active = True
        messagebox.showinfo(title="Imported Targets", message= "There were " + str(total_read) + " targets read and " + str(total_imported) + " targets imported")
    
    def process_csv_file(self, lines, total_read, total_imported):
        del lines[0]
        prev_target = ""
        text_variable = []
        for line in lines:
            parts = line.split(",")
            target = "'" + parts[0][:-2].replace("'","`") + "'"
            if target != prev_target and len(text_variable) > 0:
                total_read, total_imported = self.insert_imported_target(text_variable, total_read, total_imported)
                text_variable = []
                prev_target = ""
            if len(prev_target) == 0: 
                found = False
                for t in range(len(self.text_variables)):
                    if target == self.text_variables[t][0][0]:
                        found = True
                        break
                if found: 
                    status = messagebox.askyesnocancel("Target already exists", str(self.text_variables[t]) + "\nDo you want to Add(Yes), Overwrite(No), Skip(Cancel)")
                    if status == False: del self.text_variables[t]
                    if status == None: continue
            text_variable = self.create_ipac_entries(text_variable, target, parts, len(prev_target) == 0)
            prev_target = target
        return self.insert_imported_target(text_variable, total_read, total_imported)
        
    def insert_imported_target(self, text_variable, total_read, total_imported):
        if len(text_variable) > 0:
            text_variable.append(("",))
            total_read += 1
            if len([item for item in text_variable if 'X' in item[0]]) == 0:
                total_imported += 1
                self.text_variables.append(text_variable)
            else:
                print("Target not processed missing data: " + str(text_variable))
        return total_read, total_imported
    
    def create_ipac_entries(self, text_variable, target, parts, do_steller_entries):
        if do_steller_entries:
            text_variable.append((target,))
            text_variable.append((parts[24],) if len(parts[24]) else ('X',))
            text_variable.append((str("{:.3f}".format((float(parts[25]) + abs(float(parts[26])))/2)),) if len(parts[25]) and len(parts[26]) else ('X',))
            text_variable.append((parts[28],) if len(parts[28]) else ('X',))
            text_variable.append((str("{:.3f}".format((float(parts[29]) + abs(float(parts[30])))/2)),) if len(parts[29]) and len(parts[30]) else ('X',))
            text_variable.append((parts[20],) if len(parts[20]) else ('X',))
        text_variable.append((parts[16] if len(parts[16]) else "'?'",))
        text_variable.append(("(" + (parts[7] if len(parts[7]) else 'X') + ", 'Radius')" if parts[10] == 0 else "(" + (parts[11] if len(parts[11]) else 'X') + ", 'Mass')",))
        text_variable.append((parts[3],) if len(parts[3]) else ("X",))
        text_variable.append(("'" + parts[0][-1] + "'",))
        return text_variable

    def load_file(self, filename = None, overwrite = False):
        filename = filedialog.askopenfiles(title="Open " + ("Overwrite" if overwrite else "Append"), initialfile = self.loaded_filename, mode = 'r') if filename == None else [open(filename, 'r')]
        if filename == None or len(filename) == 0: return
        self.validate_active = False
        try:
            len(self.text_variables)
        except:
            self.text_variables = []
        for f in filename:
            self.loaded_filename = f.name
            if ".csv" in f.name:
                self.load_csv_file(f, overwrite)
            else:
                self.load_py_file(f, overwrite)
        self.undo_list.append(copy.deepcopy(self.text_variables))
        self.entry_boxes[0][0][0].focus_set()
        self.display_rows()
        self.validate_active = True
            
    def load_py_file(self, filename, overwrite):
        self.lines = filename.readlines()
        filename.close()
        start = False
        done = False
        for line in self.lines:
            if start and line.strip() == '}':
                done = True
            if start and not done:
                (line, _, comments) = line.partition('#')
                comments = comments.replace(",",";")
                line =  re.sub('\s+',' ', line).strip()
                if len(line) > 0:
                    try:
                        text_variable = []
                        (target, _, value) = line.partition(':')
                        text_variable.append((target,))
                        parts = value[3:-2].replace("[","").split(']')
                        parts1 = parts[0].split(",")
                        for v in parts1:
                            text_variable.append((v.strip(),))
                        for p in range(1, len(parts)):
                            if len(parts[p]) > 0:
                                count=0
                                for p1 in re.sub(r'\([^)]+\)', lambda x: x.group().replace(', ', '^ '), parts[p].lstrip(", ")).split(","): 
                                    count+=1
                                    text_variable.append((p1.replace('"',"").strip().replace('^ ',', '),))
                        text_variable.append((' #' + comments.replace('"', '').rstrip() if len(comments) > 0 else "",))
                        if overwrite:
                            for i in self.text_variables:
                                if i[0] == text_variable[0]:
                                    self.text_variables[i] = text_variable
                        else:
                            self.text_variables.append(text_variable)
                    except:
                        print(line)
                        print(traceback.format_exc())
            elif line.find('self.targets = {') != -1:
                start = True

    def load_csv_file(self, csv_filename, overwrite):
        csv_lines = csv_filename.readlines()
        csv_filename.close()
        csv_lines[0]=csv_lines[0].replace("\n","")
        if csv_lines[0].count(',') != len(self.headers_info)-1 and len(self.text_variables) > 0:
            parts = csv_lines[0].split(",")
            for target_field in range(len(self.headers_info)):
                if self.headers_info[target_field][3] == parts[0]:
                    break
            for name_field in range(len(self.headers_info)):
                if self.headers_info[name_field][3] == parts[1]:
                    break
            for field in range(len(self.headers_info)):
                if self.headers_info[field][3] == parts[2]:
                    break
            for line in csv_lines[1:]:
                line=line.replace("\n","")
                parts=line.split(",")
                found = False
                for tv in self.text_variables:
                    inx1 = field
                    for inx in range(name_field, len(tv), self.planet_fields):
                        if parts[0] + parts[1] == tv[target_field][0] + tv[inx][0]:
                            tv[inx1]=(parts[2].replace("'",""),)
                            found = True
                            break
                        inx1 += self.planet_fields
                if not found:
                    messagebox.showerror(title="Target not found", message= "Target " + parts[0].replace('"','') + " Planet " + parts[1].replace('"','') + " not in editor no data updated.")
        else:
            del csv_lines[0]
            text_variable = ""
            text_variables = []
            for line in csv_lines:
                if line[0] != ",":
                    if len(text_variable) > 0:
                        text_variables.append(text_variable)
                        text_variable = ""
                else:
                    line = line[self.planet_offset:]
                text_variable = text_variable + line.strip().replace('"','')
            if len(text_variable) > 0: text_variables.append(text_variable)
            for line in text_variables:
                text_variable = []
                parts = re.sub(r'\([^)]+\)', lambda x: x.group().replace(',', '^ '), line).split(",")
                for x in parts:       
                    text_variable.append((re.sub('\s+',' ', x.replace('^ ',', ').strip()),))
                if overwrite:
                    for i in range(len(self.text_variables)):
                        if self.text_variables[i][0][0] == text_variable[0][0]:
                            self.text_variables[i] = text_variable
                else:
                    self.text_variables.append(text_variable)
        
    def sort_column(self, event):
        for l in range(len(self.label_list)):
            if event.widget == self.label_list[l]:
                self.text_variables = sorted(self.text_variables, key=lambda x: tuple(x[l][0].lower()))
                break
        if self.reverse:
            self.text_variables.reverse()
        self.reverse = not self.reverse
        self.display_rows()
        
    def display_rows(self, y1 = None, x1 = None, focus_widget = None, bg = "SlateGray1", changes = True):
        self.validate_active = False
        self.title = "Targets Dictionary --> " + (self.loaded_filename if self.loaded_filename != None else "") + "    [" + (str(len(self.text_variables)) if self.text_variables != None else "0") + "]"
        self.root.title(self.title)
        if changes: self.check_for_changes()
        if y1 == None:
            focus_widget = self.root.focus_get()
            self.entry_boxes[-1][-1][0].focus_set()
        for row in range(len(self.entry_boxes)):
            for col in range(self.headers_length):
                self.entry_boxes[row][col][0].delete(0, 'end')
                self.entry_boxes[row][col][0].config(state='readonly')
                self.entry_boxes[row][col][0].config(readonlybackground ='light grey')
        row = 0
        try:
            self.bottom_row = self.top_row
            for y in range(len(self.entry_boxes)):
                bg = 'SlateGray1' if bg == 'white' else "white" if bg == "SlateGray1" else "SlateGray1"
                col = 0
                self.bottom_row += 1
                for x in range(len(self.text_variables[self.top_row + y])): 
                    if col == self.headers_length -1 and  x < len(self.text_variables[self.top_row + y]) -1:
                        col = self.planet_offset
                        row += 1
                    self.entry_boxes[row][col][0].config(state='normal', bg = bg if focus_widget == None or self.entry_boxes[row][col][0] != focus_widget else "yellow")
                    self.entry_boxes[row][col][0].insert(0, self.text_variables[self.top_row + y][x][0])
                    self.entry_boxes[row][col] = (self.entry_boxes[row][col][0], self.top_row + y, x)
                    if y1 != None and self.entry_boxes[row][col][1] == y1 and self.entry_boxes[row][col][2] == x1: 
                        focus_widget = self.entry_boxes[row][col][0]
                        threading.Thread(target=self.highlight_find, args=(row,col)).start()
                    col += 1
                row += 1
        except: pass
        if focus_widget != None: focus_widget.focus_set()
        self.validate_active = True
    
    def highlight_find(self, row, col):
        bg_color = self.entry_boxes[row][col][0].cget("background")
        self.entry_boxes[row][col][0].config(background ='yellow')
        time.sleep(1)
        self.entry_boxes[row][col][0].config(background = bg_color)
        
    def find_widget(self, widget):
        for e in range(len(self.entry_boxes)):
            for w in range(len(self.entry_boxes[e])):
                if widget == self.entry_boxes[e][w][0]:
                    return e, w
    
    def move_display(self, move_type):
        if self.text_variables == None: return
        if move_type == "scroll up":
            e, w = self.find_widget(self.root.focus_get())
            if e > 0:
                self.entry_boxes[e-1][w][0].focus_set()
            else:
                self.top_row = self.top_row - 1 if self.top_row > 0 else 0
        elif move_type == "scroll dn":
            e, w = self.find_widget(self.root.focus_get())
            if e < len(self.entry_boxes)-1:
                self.entry_boxes[e+1][w][0].focus_set()
            else:
                self.top_row = self.top_row + 1 if self.bottom_row < len(self.text_variables) +1 else self.top_row
        elif move_type == "page up":
            self.top_row = self.top_row - self.max_rows if self.top_row - self.max_rows > 0 else 0
        elif move_type == "page dn":
            self.top_row = (self.top_row + self.max_rows) if self.bottom_row + self.max_rows < len(self.text_variables) else self.top_row + (len(self.text_variables) - self.bottom_row +1)
        elif move_type == "home":
            self.top_row = 0
        elif move_type == "end":
            self.top_row = len(self.text_variables) -10
        self.display_rows()
    
    def check_for_changes(self):
        changed_data = False
        if self.text_variables != None and len(self.undo_list) > 0:
            if self.text_variables != self.undo_list[0]:
                self.undo_list.append(copy.deepcopy(self.text_variables))
                changed_data = True
            try: self.undo_button.configure(state ='normal' if len(self.undo_list) > 1 else 'disabled')
            except: pass
            self.title = self.title.replace(" *", "") + (" *" if changed_data else "")
            self.root.title(self.title)
        
    def focus(self, event):
        if self.previous_widget != None:
            self.previous_widget.config(background = self.previous_color[0])
            self.previous_widget.config(readonlybackground = self.previous_color[1])
        e, w = self.find_widget(event.widget)
        self.find_loc = self.entry_boxes[e][w][1], self.entry_boxes[e][w][2]
        for row in range(len(self.entry_boxes)):
            for col in range(self.headers_length):
                if self.entry_boxes[row][col][0].cget('readonlybackground') == "yellow":
                    self.entry_boxes[row][col][0].config(readonlybackground ='light grey')
        self.previous_color = (event.widget.cget('background'),  event.widget.cget('readonlybackground'))
        event.widget.config(background = 'yellow')
        event.widget.config(readonlybackground = 'yellow')
        self.previous_widget = event.widget
                
    def do_undo(self):
        if len(self.undo_list) > 1:
            self.text_variables = copy.deepcopy(self.undo_list[-2])
            del self.undo_list[-1]
            self.display_rows(changes = False)
        self.title = self.title.replace(" *", "") + (" *" if self.text_variables != self.undo_list[0] else "")
        self.root.title(self.title)
        
    def on_closing(self):
        if self.close_program: return
        self.close_program = True
        self.remove_find_popup()
        if "*" in self.title and messagebox.askyesno("Save Targets Data", "Do you want to save your changes?"):
            self.save_file()
        self.root.destroy()
        self.rootMain.destroy()
                
    def save_file(self):
        filename = filedialog.asksaveasfile(initialfile= self.loaded_filename, mode = 'w')
        if (filename != None):
            if ('.csv' in filename.name):
                for i in range(len(self.headers_info)):
                    filename.write(self.headers_info[i][3] + ("," if i < len(self.headers_info)-1 else ""))
                filename.write("\n")
                for line in self.text_variables:
                    try:
                        csv_line = ""
                        l_cnt = 0
                        p_cnt = 0
                        for l in range(len(line)):
                            value = "".join(line[l])
                            if (len(value) > 0):
                                if value[0] == "(":  value = '"' + value
                                if value[-1] == ")": value += '"'
                                csv_line += value + ","
                                l_cnt += 1
                                if l_cnt > self.stellar_fields + 1:
                                    p_cnt += 1
                                    if p_cnt == self.planet_fields and l < len(line) -2:
                                        csv_line += "\n,,,,,,"
                                        p_cnt = 0
                        filename.write(csv_line[:-1] + '\n')
                    except:
                        print(line)
                        print(traceback.format_exc())
            else:
                found = False
                headers_found = False
                try:
                    len(self.lines)
                except:
                    fname = filedialog.askopenfile(title='Please Choose A Targets File to Act as a Template for the Python Code.',mode = 'r')
                    if fname == None: return
                    self.lines = fname.readlines()
                    fname.close()
                for line in self.lines:
                    if line.find("# COLUMN HEADINGS:") != -1:
                        filename.write("        # COLUMN HEADINGS: ")
                        for i in range(len(self.headers_info)):
                            filename.write(self.headers_info[i][3] + (", " if i < len(self.headers_info) -1 else ""))
                        filename.write("\n")
                        headers_found = True
                        continue
                    if line.find('self.targets = {') == -1 and not found:
                        filename.write(line)
                    elif line.find('self.targets = {') != -1:
                        if not headers_found:
                            filename.write("        # COLUMN HEADINGS: ")
                            for i in range(len(self.headers_info)):
                                filename.write(self.headers_info[i][3] + (", " if i < len(self.headers_info) -1 else ""))
                            filename.write("\n")
                        filename.write(line)
                        found = True
                        for csv_line in self.text_variables:
                            if len(csv_line) > 0:
                                comment = ""
                                data = ""
                                index = 0
                                for i in range(len(csv_line)):
                                    d = csv_line[i][0].strip()
                                    if d.find("#") != -1:
                                        comment = d
                                    elif len(d) > 0:
                                        data = data + self.headers_info[index][0] + d + self.headers_info[index][1]
                                    index += 1
                                    if index == self.headers_length -1:
                                        index = self.planet_offset
                                filename.write(' '*8 + data[:-1] + ")," + ("" if len(comment) == 0 else " " + comment) + "\n")
                    elif line.strip() == '}':
                        found = False
                        filename.write(line)
                    elif line.strip()[0] == "#":
                        filename.write(line)
            filename.close()
            self.title = self.title.replace(" *", "")
            self.root.title(self.title)
        
    def modify_entry(self, widget, value, text):
        if self.validate_active:
            for e in self.entry_boxes:
                for w in e:
                    if widget == str(w[0]):
                        self.text_variables[w[1]][w[2]] = (value,)
            self.check_for_changes()
        return True

    def delete_target(self, widget, clear = False):
        e, w = self.find_widget(widget)
        if clear:
            self.text_variables[self.entry_boxes[e][w][1]] = [(' ',) for _ in range(len(self.text_variables[self.entry_boxes[e][w][1]]))]
        else:
            del self.text_variables[self.entry_boxes[e][w][1]]
        self.display_rows()
        
    def insert_target(self, widget):
        e, w = self.find_widget(widget)
        self.text_variables.insert(self.entry_boxes[e][w][1], [(' ',) for _ in range(self.headers_length)])
        self.display_rows()
        for col in range(self.headers_length):
            self.entry_boxes[e][col][0].config(state='normal')
                               
    def clicker_target(self, e):
        rmenu = tk.Menu(None, tearoff = 0, takefocus = False) 
        rmenu.add_command(label = 'Insert Target', command = (lambda: self.insert_target(e.widget))) 
        rmenu.add_command(label = 'Delete Target', command = (lambda: self.delete_target(e.widget))) 
        rmenu.add_command(label = 'Clear target', command = (lambda: self.delete_target(e.widget, True)))
        rmenu.tk_popup(e.x_root + 40, e.y_root + 10, entry = "0")
    
    def delete_planet(self, widget, clear = False):
        e, w = self.find_widget(widget)
        if len(self.text_variables[self.entry_boxes[e][w][1]]) > self.headers_length:
            for j in range(self.planet_fields):
                if clear:
                    self.text_variables[self.entry_boxes[e][w][1]][self.entry_boxes[e][w][2]+j] = (' ',)
                else:
                    del self.text_variables[self.entry_boxes[e][w][1]][self.entry_boxes[e][w][2]]
            self.display_rows()
         
    def insert_planet(self, widget):
        e, w = self.find_widget(widget)
        if len(self.text_variables[self.entry_boxes[e][w][1]]) > self.headers_length -1:
            for _ in range(self.planet_fields):
                self.text_variables[self.entry_boxes[e][w][1]].insert(self.entry_boxes[e][w][2], (' ',))
            self.display_rows()
            for col in range(self.planet_fields):
                self.entry_boxes[e][w + col][0].config(state='normal')
        
    def clicker_planet(self, e):
        rmenu = tk.Menu(None, tearoff = 0, takefocus = False) 
        rmenu.add_command(label = 'Insert Planet', command = (lambda: self.insert_planet(e.widget))) 
        rmenu.add_command(label = 'Delete Planet', command = (lambda: self.delete_planet(e.widget)))
        rmenu.add_command(label = 'Clear Planet', command = (lambda: self.delete_planet(e.widget, True))) 
        rmenu.tk_popup(e.x_root + 40, e.y_root + 10, entry = "0")
    
    def create_help_popup(self):
        if not hasattr(self, "help_top") or not self.help_top.winfo_exists():
            self.help_top = tk.Toplevel(self.root)
            self.help_top.attributes('-topmost', 'true')
            self.help_top.title('Help')
            tk.Label(self.help_top, text="Cursor Down -> <Down Arrow>").grid(row=0, column=0, sticky='w', padx =5, pady = 0)
            tk.Label(self.help_top, text="Cursor Up -> <Up Arrow>").grid(row=0, column=1, sticky='w', padx =5, pady = 0)
            tk.Label(self.help_top, text="Page Down -> <PgDn>").grid(row=0, column=2, sticky='w', padx =5, pady = 0)
            tk.Label(self.help_top, text="Page Up -> <PgUp>").grid(row=1, column=0, sticky='w', padx =5, pady = 0)
            tk.Label(self.help_top, text="Cursor End -> <End>").grid(row=1, column=1, sticky='w', padx =5, pady = 0)
            tk.Label(self.help_top, text="Cursor Home -> <Home>").grid(row=1, column=2, sticky='w', padx =5, pady = 0)
            tk.Label(self.help_top, text="Cursor Right -> <Tab>").grid(row=2, column=0, sticky='w', padx =5, pady = 0)
            tk.Label(self.help_top, text="Cursor Left -> <Shift Tab>").grid(row=2, column=1, sticky='w', padx =5, pady = 0)
            tk.Label(self.help_top, text="Find -> <Control-f>").grid(row=2, column=2, sticky='w', padx =5, pady = 0)
            tk.Label(self.help_top, text="Load targets Append-> <Control-l>").grid(row=3, column=0, sticky='w', padx =5, pady = 0)
            tk.Label(self.help_top, text="Undo -> <Control-z>").grid(row=3, column=1, sticky='w', padx =5, pady = 0)
            tk.Label(self.help_top, text="Help -> <F1>").grid(row=3, column=2, sticky='w', padx =5, pady = 0)
            tk.Label(self.help_top, text="Remove Popups -> <Esc>").grid(row=4, column=0, sticky='w', padx =5, pady = 0)
            tk.Label(self.help_top, text="Display Popups -> <Mouse Button-3>").grid(row=4, column=1, sticky='w', padx =5, pady = 0)
            tk.Label(self.help_top, text="Save Targets -> <Control-s>").grid(row=4, column=2, sticky='w', padx =5, pady = 0)
            tk.Label(self.help_top, text="Import Targets -> <Control-i>").grid(row=5, column=0, sticky='w', padx =5, pady = 0)
            tk.Label(self.help_top, text="Delete All Targets -> <Control-Delete>").grid(row=5, column=0, sticky='w', padx =5, pady = 0)
            tk.Label(self.help_top, text="Modify Header -> <Control-m>").grid(row=5, column=1, sticky='w', padx =5, pady = 0)
            tk.Label(self.help_top, text="Load targets Overwrite-> <Control-L>").grid(row=5, column=2, sticky='w', padx =5, pady = 0)
            self.close_button = tk.Button(self.help_top, text="Close", command=lambda event: self.remove_find_popup())
            self.close_button.grid(row=6, column=1, padx = 5, pady = 0)
            self.close_button.bind("<Return>", lambda event: self.remove_find_popup())
            self.close_button.bind("<Escape>", lambda event: self.remove_find_popup())
            self.close_button.focus_set()
        
    def create_find_popup(self):
        if not hasattr(self, "find_top") or not self.find_top.winfo_exists():
            self.find_top = tk.Toplevel(self.root)
            self.find_top.attributes('-topmost', 'true')
            self.find_top.title('Find/Replace')
            tk.Label(self.find_top, text="Find:").grid(row=0, column=0, sticky='w', padx =10, pady = 10)
            self.find_entry = ttk.Combobox(self.find_top, width = 15, justify = "left", validate = "key", validatecommand = self.vcmd2)
            self.find_entry["values"] = self.find_values_list
            self.find_entry.configure(state = 'normal')
            self.find_entry.bind("<Return>", self.find)
            self.find_entry.bind("<Escape>", lambda event: self.remove_find_popup())
            self.find_entry.bind('<<ComboboxSelected>>', self.enable_find_buttons)
            self.find_entry.grid(row=0, column=1, sticky='w', padx = 10, pady = 10, columnspan = 3)
            tk.Label(self.find_top, text="Replace:").grid(row=1, column=0, sticky='w', padx = 10, pady = 10)
            self.replace_entry = ttk.Combobox(self.find_top, width = 15, justify = "left", validate = "key", validatecommand = self.vcmd1)
            self.replace_entry["values"] = self.replace_values_list
            self.replace_entry.configure(state = 'normal')
            self.find_entry.bind("<Return>", self.replace)
            self.replace_entry.bind("<Escape>", lambda event: self.remove_find_popup())
            self.replace_entry.bind('<<ComboboxSelected>>', self.enable_replace_buttons)
            self.replace_entry.grid(row=1, column=1, sticky='w', padx =10, pady = 10, columnspan = 3)
            self.find_button = tk.Button(self.find_top, text='Find', command=self.find, state='disable')
            self.undo_button = tk.Button(self.find_top, text='Undo', command=self.undo, state='disable')
            self.undo_button.grid(row=2, column=2, sticky='w', padx = 10, pady = 0)
            self.find_button.grid(row=2, column=0, sticky='e', padx = 20, pady = 10)
            self.replace_button = tk.Button(self.find_top, text='Replace', command=self.replace, state='disable')
            self.replace_button.grid(row=3, column=0, sticky='w', padx = 10, pady = 10)
            self.find_replace_button = tk.Button(self.find_top, text='Replace/Find', command=self.find_and_replace, state='disable')
            self.find_replace_button.grid(row=2, column=1, sticky='w', padx = 10, pady = 10)
            self.replace_all_button = tk.Button(self.find_top, text='Replace All', command=(lambda: self.replace(None, True, True)), state='disable')
            self.replace_all_button.grid(row=3, column=1, sticky='w', padx = 10, pady = 0)
            tk.Button(self.find_top, text='Close', command=self.remove_find_popup).grid(row=3, column=2, sticky='w', padx = 10, pady = 0)
            self.status_label = tk.Label(self.find_top, text=" ")
            self.status_label.grid(row=4, column=0, sticky='w', padx =0, pady =0, columnspan = 3)
            self.find_top.geometry("+%d+%d" % (self.root.winfo_x() + 600, self.root.winfo_y() + 200))
            self.find_entry.focus_set()
    
    def undo(self):
        self.status_label.config(text="")
        self.do_undo()
        
    def remove_find_popup(self):
        try: self.find_top.destroy()
        except: pass
        try: self.help_top.destroy()
        except: pass
        try: self.modify_header_top.destroy()
        except: pass
    
    def enable_find_buttons(self, value):
        if not isinstance(value, str): value = self.find_entry.get()
        self.find_button.configure(state ='normal' if len(value) > 0 else 'disabled')
        return True

    def enable_replace_buttons(self, value = None):
        if not isinstance(value, str): value = self.replace_entry.get()
        self.replace_button.configure(state ='normal' if len(value) > 0 else 'disabled')
        self.find_replace_button.configure(state ='normal' if len(value) > 0 else 'disabled')
        self.replace_all_button.configure(state ='normal'  if len(value) > 0 else 'disabled')
        return True
            
    def find_and_replace(self):
        self.replace()
        self.find()
        
    def replace(self, event = None, clear_status = True, replace_all = False):
        if clear_status: self.status_label.config(text="")
        value = self.find_entry.get()
        value1 = self.replace_entry.get()
        self.find_replace(value, value1, replace_all)
         
    def find(self, event = None, clear_status = True):
        if clear_status: self.status_label.config(text="")
        value = self.find_entry.get().lower()
        self.find_replace(value)
    
    def find_replace(self, value, value1 = None, replace_all = False):
        replace_count = 0
        if value not in self.find_values_list:
            self.find_values_list.append(value)
            self.find_entry["values"] = self.find_values_list
        if value1 != None and value1 not in self.replace_values_list:
            self.replace_values_list.append(value1)
            self.replace_entry["values"] = self.replace_values_list
        prev_find_loc = self.find_loc
        x_offset = self.find_loc[1]
        for y in range(self.find_loc[0], len(self.text_variables)):
            for x in range(x_offset, len(self.text_variables[y])): 
                if (self.text_variables[y][x][0].lower() == value) if replace_all else value in self.text_variables[y][x][0].lower() and (self.find_loc != (y,x) or value1 != None) or self.text_variables[y][x][0].lower() == value and self.find_loc == (0, 0):
                    if value1 != None:
                        self.text_variables[y][x] = (value1,)
                        replace_count += 1
                    self.top_row = y
                    if not replace_all:
                        self.display_rows(y,x)
                        return
            x_offset = 0
        if replace_all:
            self.display_rows()
            self.status_label.config(text="Replaced " + str(replace_count) + " Occurrences")
        elif self.find_loc == prev_find_loc:
            if self.find_loc == (0, 0):
                self.status_label.config(text="Data not found!!")
                return
            self.find_loc = (0, 0)
            return self.find_replace(value, value1, False)
        
if __name__ == "__main__":
    if len(sys.argv) != 3: 
        root=tk.Tk()
        root.withdraw()
    #headerxx=[('', ': ', 'T', 'TARGET', '15'), ('([', ', ', 'S', 'RADIUS', '8'), ('', ', ', 'S', 'RADIUS~', '9'), ('', ', ', 'S', 'MASS', '7'), ('', ', ', 'S', 'MASS~', '7'), ('', '],', 'S', 'TEMP', '7'), (' [', ', ', 'P', 'PERIOD', '20'), ('', ', ', 'P', 'RADIUS/MASS', '30'), ('', ', ', 'P', 'INCLINATION', '28'), ('', ', ', 'P', 'ECCENTRICITY', '14'), ('', '],', 'P', 'NAME', '6'), ('', '', 'C', 'COMMENTS', '35')]
    headerxx=[('', ': ', 'T', 'TARGET', '15'), ('([', ', ', 'S', 'RADIUS', '8'), ('', ', ', 'S', 'RADIUS~', '9'), ('', ', ', 'S', 'MASS', '7'), ('', ', ', 'S', 'MASS~', '7'), ('', '],', 'S', 'TEMP', '7'), (' [', ', ', 'P', 'PERIOD', '20'), ('', ', ', 'P', 'RADIUS/MASS', '30'), ('', ', ', 'P', 'INCLINATION', '28'), ('', ', ', 'P', 'ECCENTRICITY', '14'), ('', '],', 'P', 'NAME', '6'), ('', '', 'C', 'COMMENTS', '35')]
    dynamite_targets_editor(sys.argv[0] if len(sys.argv) == 3 else root, sys.argv[1] if len(sys.argv) >= 2 else None, headerxx, None)


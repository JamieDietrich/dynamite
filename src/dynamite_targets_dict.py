### DYNAmical Multi-planet Injection TEster (DYNAMITE) ###
### Targets Dictionary-CSV Conversion ###
### Jeremy Dietrich ###
### jdietrich1@email.arizona.edu ###
### 2020 August 12 ###
### Version 1.3 ###
### Dietrich & Apai (2020), Astronomical Journal ###
### https://doi.org/10.3847/1538-3881/aba61d ###

import re
import sys
import traceback

csv_headers = "TARGET,STELLAR RADIUS,STELLAR RADIUS UNCERTAINTY,STELLAR MASS,STELLAR MASS UNCERTAINTY,STELLAR TEMPERATURE,PLANET INCLINATION,PLANET RADIUS/MASS,PLANET PERIOD,PLANET NAME,COMMENTS\n"

class dynamite_targets_dict:
    def __init__(self, infile, outfile):
        """ Check input file types"""
        if infile.find('.csv') != -1 and outfile.find('.py') != -1:
            print("Running CSV2PY for", infile, "to", outfile + ".new")
            self.csv_to_targets(infile, outfile)
            print("Complete CSV2PY for", infile, "to", outfile + ".new")
        elif infile.find('.py') != -1 and outfile.find('.csv') != -1:
            print("Running PY2CSV for", infile, "to", outfile)
            self.targets_to_csv(infile, outfile)
            print("Complete PY2CSV for", infile, "to", outfile)
        else:
            print('Input files must consist of a python targets file and a tragets csv file!')

    def targets_to_csv(self, py_filename, csv_filename):
        """Targets dict to csv file"""
        with open(py_filename, 'r') as f:
            lines = f.readlines()
        start = False
        done = False
        with open(csv_filename, mode='w', newline='') as csv_file:
            csv_file.write(str(csv_headers))
            for line in lines:
                if start and line.strip() == '}':
                    done = True
                if start and not done:
                    (line, _, comments) = line.partition('#')
                    line = line.strip()
                    if len(line) == 0:
                        continue
                    try:
                        (csv_line, _, value) = re.sub('\s+',' ', line).partition(':')
                        parts = value[3:-2].replace("[","").split(']')
                        csv_line += "," + parts[0].replace(" ","") + parts[1] + ",\n"
                        for p in range(2, len(parts)):
                            if len(parts[p]) > 0:
                                csv_line +=  ",,,,," + parts[p] + ",\n" 
                                csv_line = csv_line.replace(' (', '"(').replace('),', ')",').replace("',\"", "',")
                        csv_file.write(" " + csv_line[:-1] + ('"  #' + comments.replace('"', '').rstrip() + '"\n' if len(comments)> 0 else "\n"))
                    except:
                        print(line)
                        print(traceback.format_exc())
                elif line.find('self.targets = {') != -1:
                    start = True

    def csv_to_targets(self, csv_filename, py_filename):
        """CSV to Targets"""
        with open(py_filename, 'r') as f:
            lines = f.readlines()
        with open(py_filename + ".new", mode='w') as pyfile:
            found = False
            for line in lines:
                if line.find('self.targets = {') == -1 and not found:
                    pyfile.write(line)
                elif line.find('self.targets = {') != -1:
                    pyfile.write(line)
                    found = True
                    with open(csv_filename, newline = '') as csvfile:
                        csv_lines = csvfile.readlines()
                    row = 1
                    while row < len(csv_lines):
                        csv_lines[row] = re.sub('\s+',' ', csv_lines[row])
                        parts = csv_lines[row][:-1].split(',', 6)
                        current_row = "'" + parts[0] + ": ([" +', '.join([parts[1].lstrip().replace("'",""),[', '.join(parts[2:-1])][0]]) + "], [" + parts[-1].lstrip()         
                        current_row = current_row.replace('"', '')[:-1].replace("),(", "), (") + "],"
                        row += 1
                        while row < len(csv_lines) and csv_lines[row][0] == ",":
                            value, _, comments = csv_lines[row].partition("#")
                            current_row += ' [' + value.rstrip().replace('"', '')[6:-1].lstrip().replace("),(", "), (") + '],'
                            row += 1
                        current_row = current_row.replace('"','').replace("', ", ", ").replace(", ' ",", ").replace("', (", ", (").replace("',(", "', (").replace(",(", ", (")
                        pyfile.write(' '*7 + current_row[1: -1] + ")," + (" #" + comments.replace('"','').rstrip() if len(comments) > 0 else "") + '\n')
                elif line.strip() == '}':
                    found = False
                    pyfile.write(line)
                elif line.strip()[0] == "#":
                    pyfile.write(line)

if __name__ == '__main__':
    if len(sys.argv) == 3:
        dynamite_targets_dict(sys.argv[1], sys.argv[2])
    else:
        print("Not enough input files")
        print('Usage: One Source and One Destination file ie. dynamite_targets.py  targets.csv or targets.csv  dynamite_targets.py')

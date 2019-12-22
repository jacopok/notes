import os
import re

align_start = "\\begin{align}"
align_end = "\\end{align}"
subeq_start = "\\begin{subequations}"
subeq_end = "\\end{subequations}"
align_newline = "\\\\"

folder_names = [
"advanced_astrophysics",
"astrophysics_lab",
"general_relativity",
"numerical_methods",
"astrophysics_cosmology"
]
filenames_regex = "\\d{2}\\w{3}\\.tex"

def make_subeq_align(filename):
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    for i, line in enumerate(lines):
        if (align_start in line
            and
            subeq_start not in lines[i - 1]):
            
            newline = False
            for pnumber, pline in enumerate(lines[i:]):
                if align_newline in pline:
                    newline = True
                if align_end in pline:
                    j = pnumber + i
                    break
            if (not newline):
                pass
            else:
                lines[i] = subeq_start + '\n' + lines[i]
                lines[j] = lines[j] + subeq_end + '\n'
    
    with open(filename, 'w') as f:
        f.writelines(lines)

for folder in folder_names:
    for fname in os.listdir(folder):
        if (re.match(filenames_regex, fname)):
            make_subeq_align(os.path.join(folder, fname))
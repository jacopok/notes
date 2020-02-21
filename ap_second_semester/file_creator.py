import os
import re
import datetime

subfile_heading = """\\documentclass[main.tex]{subfiles}
\\begin{document}

\\end{document}"""

end_doc = "\\end{document}\n"
subfile_include = lambda x : "\\subfile{" + x + "}\n"

folder_names = {
    # "aa": "advanced_astrophysics",
    # "al": "astrophysics_lab",
    # "gr": "general_relativity",
    # "nm": "numerical_methods",
    # "ac": "astrophysics_cosmology"
}

schedule = {
    0: [],
    1: [],
    2: ["ac"],
    3: ["gr", "ac"],
    4: ["gr"]
}

def next_monday():
    """
    Returns a datetime object for next monday if it is
    Wednesday or later,
    for today or yesterday if it is monday or tuesday
    """
    now = datetime.date.today()
    weekday = now.weekday()
    if (weekday == 0):
        return now
    elif (weekday == 1):
        return now - datetime.timedelta(1)
    elif (weekday > 1):
        return now + datetime.timedelta(7 - weekday)

def date_filename(date):
    return (date.strftime("%b%d").lower())
    
def create_file(folder, filename):
    main_path = os.path.join(folder, "main.tex")
    file_path = os.path.join(folder, filename + ".tex")

    with open(file_path, "x") as f:
        f.write(subfile_heading)
        print("Creating file at " + file_path)

    with open(main_path, "r") as f:
        buf = f.readlines()

    with open(main_path, "w") as f:
        for line in buf:
            if line == end_doc:
                line = subfile_include(filename) + line
                print("Adding subfile to main at " + main_path)
            f.write(line)

def add_line_main(folder):
    main_path = os.path.join(folder, "main.tex")

    with open(main_path, "r") as f:
        buf = f.readlines()

    with open(main_path, "w") as f:
        for line in buf:
            if line == end_doc:
                line = "\n" + line
                print("Adding white line to main at " + main_path)
            f.write(line)

for i in schedule:
    day = next_monday() + datetime.timedelta(i)
    filename = date_filename(day)

    lessons = schedule[i]
    
    for lesson in lessons:
        folder = folder_names[lesson]
        create_file(folder, filename)
    
for folder in folder_names:
    add_line_main(folder_names[folder])

import json
import os
import pywavefront  
import numpy as np  
from plyfile import PlyData, PlyElement  

def read_lines_from_file(file_path):
    lines = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('-'):
                continue
            floats = line.strip().split()
            if float(floats[3]) < 0.2 and float(floats[5])>0.8 :
                lines.append(floats[:3])
    return lines

for i in range(2000):
# for i in [700]:
    file_path = f'raw/{i}.txt'
    output_file = 'tmp2/output.json'

    lines = read_lines_from_file(file_path)
    output_data = [[float(x) for x in line] for line in lines]  # Convert to list of lists
    def write_to_json(data, output_file):
        with open(output_file, 'w') as file:
            json.dump(data, file, indent=4)
    write_to_json(output_data, output_file)

    os.system(f"splashsurf reconstruct tmp2/output.json -r=0.04 -c=0.6 -l=2.4 -t=0.6 --subdomain-grid=on -o tmp3/fluid_{i}.obj")
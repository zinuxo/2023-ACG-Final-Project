import json
import os
import pywavefront  
import numpy as np  
from plyfile import PlyData, PlyElement  

def read_lines_from_file(file_path):
    lines = []
    with open(file_path, 'r') as file:
        for line in file:
            floats = line.strip().split()
            lines.append(floats[:3])
    return lines


for i in range(1,19):
    file_path = f'tmp2/tmpp{i}.txt'
    output_file = f'tmp2/output.json'
    lines = read_lines_from_file(file_path)
    output_data = [[float(x) for x in line] for line in lines]  # Convert to list of lists
    def write_to_json(data, output_file):
        with open(output_file, 'w') as file:
            json.dump(data, file, indent=4)
    write_to_json(output_data, output_file)

    os.system(f"splashsurf reconstruct tmp2/output.json -r=0.02 -c=1.5 -l=2 -t=0.6 --subdomain-grid=on -o tmp2/output{i}.obj")

# input_file = 'tmp/output.obj'
# output_file = 'tmp/output.ply'

# def obj_to_ply(obj_filename, ply_filename):  
#     # 读取.obj文件  
#     obj = pywavefront.Wavefront(obj_filename)  
#     vertices = obj.vertices
#     for mesh in obj.mesh_list:
#         print(mesh.faces)
  
#     # 创建PLY数据结构  
#     # vertex_element = PlyElement.describe(vertices, 'vertex')
#     # face_element = PlyElement.describe(meshes, 'face')
#     # data = PlyData([vertex_element, face_element])

#     # # 写入PLY文件  
#     # with open(ply_filename, 'wb') as f:  
#     #     data.write(f)  

# obj_to_ply(input_file, output_file)


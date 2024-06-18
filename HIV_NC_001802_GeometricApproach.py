import os
import tkinter
import datetime
import numpy as np
import tkinter.filedialog        
import matplotlib.pyplot as plt 
from mpl_toolkits import mplot3d
from matplotlib.backends.backend_pdf import PdfPages

os.listdir('D:\OTHERS\projects\skews n zcurves')

root = tkinter.Tk()
root.withdraw()

input_file_path = tkinter.filedialog.askopenfilename(title = "Select input file", filetypes = (("Text file","*.txt"),("FASTA file", '*.fasta'),("all files","*.*")))

if input_file_path == "":
    print('No file was chosen!')
else:
    print(input_file_path)
    print(os.path.split(input_file_path))
    print(os.path.split(input_file_path)[1])    
    print(os.path.split(input_file_path)[1].split("."))    
    print(os.path.split(input_file_path)[1].split(".")[0])
    print("\n")

    output_file_name = os.path.split(input_file_path)[1].split(".")[0] + "_RNA" + ".txt"
    output_file_path = os.path.split(input_file_path)[0] + "/" + output_file_name
    
    #output_file_path = tkinter.filedialog.asksaveasfilename(title = "Save output file", initialfile = output_file_name, filetypes = (("Text file","*.txt"),("FASTA file", '*.fasta'),("all files","*.*")))

    print("Input file: ", input_file_path)
    print("Output file: ", output_file_path)
    
temp_DNA = ''           
comp_RNA = ''           
output_file = open(output_file_path, 'w')

with open(input_file_path,'r') as input_data:
    header = input_data.readline().strip()
    output_header = header + ' RNA sequence'
    output_file.write(output_header)
    print("Input Header: ", header)
    print("Output Header: ", output_header)
    for line in input_data:
        temp_DNA = line.strip()
        comp_RNA = ""
        for base in temp_DNA:
            if base == "A":
                  comp_RNA = comp_RNA + 'U'
            elif base == "C":
                comp_RNA = comp_RNA + 'G'
            elif base == "G":
                comp_RNA = comp_RNA + 'C'
            elif base == "T":
                comp_RNA = comp_RNA + 'A'
        output_file.write("\n" + comp_RNA)
        output_file.flush()        
        print(temp_DNA)        
        print(comp_RNA)
        print("===================================================")
input_data.close()
output_file.close()

def all_variables():    
    global counter
    global A_count
    global C_count
    global G_count
    global U_count
    global gc_skew
    global au_skew
    global ca_skew
    global ga_skew
    global ua_skew
    global uc_skew   
    global ug_skew    
    global cg_skew    
    global A_percentage
    global C_percentage
    global G_percentage
    global U_percentage
    global GC_content_percentage
    global x
    global y
    global z    
    
    counter = 0    
    A_count = 0
    C_count = 0
    G_count = 0
    U_count = 0
    
    gc_skew = []
    au_skew = []    
    
    ca_skew = []
    ga_skew = []
    ua_skew = []
    uc_skew = []
    ug_skew = []
    cg_skew = []   
           
    A_percentage = []
    C_percentage = []
    G_percentage = []
    U_percentage = []
    GC_content_percentage = []
    
    x = []
    y = []
    z = []
    
    #print('done!')

def bases_skew(A, B):
    try: return (A - B) /(A + B)
    except ZeroDivisionError: return 0

print(bases_skew(0,0))

all_variables()    
with open(output_file_path,'r') as input_data:
    header = input_data.readline().strip()
    for line in input_data:
        temp_DNA = line.strip()
        for base in temp_DNA:
            counter += 1
            if base == "A":
                A_count +=1                  
            elif base == "C":
                C_count +=1
            elif base == "G":
                G_count +=1
            elif base == "U":
                U_count +=1            
            
            A_percentage.insert(counter, (A_count/counter*100))
            C_percentage.insert(counter, (C_count/counter*100))
            G_percentage.insert(counter, (G_count/counter*100))
            U_percentage.insert(counter, (U_count/counter*100))

            GC_content_percentage.insert(counter, ((G_count + C_count)/counter * 100) )

            gc_skew.insert(counter, bases_skew(G_count, C_count))
            au_skew.insert(counter, bases_skew(A_count, U_count))            

            ca_skew.insert(counter, bases_skew(C_count, A_count))
            ga_skew.insert(counter, bases_skew(G_count, A_count))
            ua_skew.insert(counter, bases_skew(U_count, A_count))
            ug_skew.insert(counter, bases_skew(U_count, G_count))
            uc_skew.insert(counter, bases_skew(U_count, C_count))
            cg_skew.insert(counter, bases_skew(C_count, G_count))            
            
            # x, y, and z are not normalized.
            x.insert(counter,((A_count + G_count)-(C_count + U_count))) # Purine vs Pyrimidine
            y.insert(counter,((A_count + C_count)-(G_count)))           # Amino vs Keto
            z.insert(counter,((A_count + U_count)-(C_count + G_count))) # Weak vs Strong Hydrogen Bonds

input_data.close()

print("Total Bases =", counter)
print("Total GC = ", C_count + G_count)
formatted_GC_ratio ="{:.2f}".format((C_count + G_count)/counter * 100)
print("GC ratio =", formatted_GC_ratio)

print("A% = {:.2f}".format(A_count/counter * 100))
print("C% = {:.2f}".format(C_count/counter * 100))
print("G% = {:.2f}".format(G_count/counter * 100))
print("U% = {:.2f}".format(U_count/counter * 100))
print(min(uc_skew))
print(max(uc_skew))

plt.figure(figsize=(9,6))
plt.plot(A_percentage, label="A%") 
plt.plot(C_percentage, label="C%") 
plt.plot(G_percentage, label="G%")
plt.plot(U_percentage, label="U%") 
plt.legend()
plt.grid()
plt.show() 
plt.close()

plt.figure(figsize=(9,6))
plt.plot(A_percentage[100:], label="A%") 
plt.plot(C_percentage[100:], label="C%") 
plt.plot(G_percentage[100:], label="G%")
plt.plot(U_percentage[100:], label="U%") 
plt.legend()
plt.grid()
plt.show() 
plt.close()

print("A% = {:.2f}".format(A_count/counter * 100))
print("\x1B[1mC% = {:.2f}".format(C_count/counter * 100))
print("\x1B[1mG% = {:.2f}".format(G_count/counter * 100))  
print("\x1b[0mU% = {:.2f}".format(U_count/counter * 100))
print("\x1b[1;35mGC% =", formatted_GC_ratio)

plt.plot(GC_content_percentage, label="GC %") 
plt.title('Total GC Percentage')
plt.tight_layout()
plt.legend()
plt.grid()
plt.show() 
plt.close()

plt.plot(GC_content_percentage[100:], label="GC%") 
plt.title('Total GC Percentage')
plt.tight_layout()
plt.legend()
plt.grid()
plt.show() 
plt.close()

plt.figure(figsize=(9,6))
plt.plot(ca_skew, label="CA Skew") 
plt.plot(ga_skew, label="GA Skew") 
plt.plot(ua_skew, label="UA Skew") 
plt.plot(ug_skew, label="UG Skew") 
plt.plot(uc_skew, label="UC Skew") 
plt.plot(cg_skew, label="CG Skew") 
plt.legend()
plt.grid()
plt.show() 
plt.close()

plt.figure(figsize=(9,6))
plt.plot(np.cumsum(ca_skew), label="Cumulative CA Skew") 
plt.plot(np.cumsum(ga_skew), label="Cumulative GA Skew") 
plt.plot(np.cumsum(ua_skew), label="Cumulative UA Skew") 
plt.plot(np.cumsum(ug_skew), label="Cumulative UG Skew") 
plt.plot(np.cumsum(uc_skew), label="Cumulative UC Skew") 
plt.plot(np.cumsum(cg_skew), label="Cumulative CG Skew") 
# We can get the virus name from the file name
virus_Name = os.path.split(input_file_path)[1].split("_")[0]
plt.title(virus_Name + " Skew Profiles")
plt.legend()
plt.grid()
plt.show() 
plt.close()

plt.figure(figsize=(8,5))
plt.plot(x, label="X: Purine (A+G) vs Pyrimidine (C+U)")
plt.plot(y, label="Y: Amino (A+C) vs Keto (G)")
plt.plot(z, label="Z: Weak (A+U) vs Strong (C+G) Hydrogen Bonds")
plt.title('Z Curve')
plt.legend()
plt.grid()
plt.show() 
plt.close()

plt.figure(figsize=(8,5))
plt.plot(x, label="X: Purine (A+G) vs Pyrimidine (C+U)") 
plt.title('Z-Curve - X axis')
plt.legend()
plt.grid()
plt.show() 
plt.close()

plt.figure(figsize=(8,5))
plt.plot(y, label="Y: Amino (A+C) vs Keto (G)", color='blue') 
plt.plot(y[:10000], label="Y: Amino (A+C) vs Keto (G) 1st part", color='orange') 
plt.plot(y[10000:20000], label="Y: Amino (A+C) vs Keto (G) 2nd part", color='red') 
plt.plot(y[20000:], label="Y: Amino (A+C) vs Keto (G) 3rd part", color='green') 
plt.title('Z-Curve - Y axis')
plt.legend()
plt.grid()
plt.show() 
plt.close()

plt.figure(figsize=(8,5))
plt.plot(z, label="Z: Weak (A+U) vs Strong (C+G) Hydrogen Bonds", color='green') 
plt.title('Z-Curve - Z axis')
plt.tight_layout()
plt.legend()
plt.grid()
plt.show() 
plt.close()

matplotlib.inline # type: ignore


fig = plt.figure(figsize=(10,8))
ax = plt.axes(projection='3d')
ax.set_xlabel('X', fontsize=20)
ax.set_ylabel('Y', fontsize=20)
ax.set_zlabel('Z', fontsize=20)
ax.scatter(x, y, z, c=z, cmap='viridis', linewidth=0.1, depthshade=True);
ax.set_title('Z Curve', fontsize=20,loc='center', pad=20)
ax.view_init(elev=30., azim=-60)
plt.tight_layout()

fig = plt.figure(figsize=(10,8))
ax = plt.axes(projection='3d')
ax.set_xlabel('X', fontsize=20)
ax.set_ylabel('Y', fontsize=20)
ax.set_zlabel('Z', fontsize=20)
ax.scatter(x, y, z, c=z, cmap='viridis', linewidth=0.1, depthshade=True);
ax.set_title('Z Curve', fontsize=20,loc='center', pad=20)
ax.view_init(elev=30., azim=-60)
plt.tight_layout()
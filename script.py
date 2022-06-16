import subprocess
import sys
import os
import shutil


def latticeParabola(x1, x2, x3, y1, y2, y3, xP, yP):
    denom = (x1 - x2)(x1 - x3)(x2 - x3)
    a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2))
    b = (x3^2 * (y1 - y2) + x2^2 * (y3 - y1) + x1^2 * (y2 - y3))
    c = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3)
    xP = -b/ (2 * a)
    yP = (c - b^2 / (4 * a)) / denom

if len(sys.argv) != 3:
    raise ValueError('Usage: arg[1] = Compound Name, arg[2] = Structure')
 
compound = sys.argv[1]
structure = sys.argv[2]

dir = "Results\\"+compound+"_"+structure
if not os.path.exists(dir):
    os.mkdir(dir)
shutil.copytree("CubicHeuslerScripts\\"+structure, dir+"\\"+structure)
os.chdir(dir)

# Lattice Optimization
# Eventually Api option
latticeEstimate = input("Enter estimate for lattice energy (https://oqmd.org):")

inputParam = str(latticeEstimate) + " " + str(latticeEstimate + 0.001) + " " + str(latticeEstimate + 0.002)
output = []
latticeEstimate = subprocess.run(["FH_lattice_optzm.bash", inputParam], stdout=subprocess.PIPE, text=True)
for line in iter(latticeEstimate.stdout.readline, b''):
    output.append(float(line))
latticeEstimate.stdout.close()
xP = 0
yP = 0
latticeParabola(latticeEstimate, latticeEstimate+ 0.001, latticeEstimate+0.002, output[0], output[1], output[2],xP, yP)

print("lattice vertex x=",xP," and y=", yP)




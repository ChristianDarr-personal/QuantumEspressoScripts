import subprocess
import sys
import os

 
if len(sys.argv) != 3:
    raise ValueError('Usage: arg[1] = Compound Name, arg[2] = Structure')
 
compound = sys.argv[1]
structure = sys.argv[2]

dir = "Results\\"+compound+'_'+structure;
if not os.path.exists(dir):
    os.mkdir(dir)






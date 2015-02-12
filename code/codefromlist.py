#------------------------------------------------------#
import sys
import numpy as N
from plotvalues import read_file_values as get_values
from fixedpoint import fpsg32 as q
import os
import os.path as path
#------------------------------------------------------#
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print "USAGE: "+path.basename(sys.argv[0])+" <raw-input-file>"
        exit(1)
    path = sys.argv[1]
    try:
        name,ext = os.path.basename(path).split(".")
    except Exception as e:
        print "\n\tERROR:\tcan only handle files with extensions!"
        print "\t\t"+str(e)
        exit(1)
    #Creating the source file#
    f = open(path)
    try:
        table = get_values(f)
    except ValueError as e:
        print "\n\tERROR:\tUnknown syntax detected!"
        print "\t\t"+str(e)
        exit(1)        
    ret = "#include \""+name+".h"+"\"\n\n"
    typ = "sg32"
    for l,key in enumerate(table):
        values = N.array(table[key])
        dim = len(values.shape)
        if dim == 1:
            ret += typ+" "+key+"["+str(len(values))+"] = "+"{\n"
            for k,i in enumerate(values):
                if k == len(values)-1:
                    ret += "\t"+str(i)+"\n"
                else:
                    ret += "\t"+str(i)+",\n"
            if l == len(table)-1:
                ret += "};"
            else:
                ret += "};\n\n"
    s = "\nValues found ["+str(len(table.keys()))+"]:\n------------------\n"
    for l,x in enumerate(table.keys()):
        if l == len(table.keys())-1:
            s += x+"."
        else:
            s += x+", "
    f.close()
    print s

    #Creating the header file#
    path,_ = path.split(".")
    path += ".c"
    f = open(path,"w+")
    f.write(ret)
    f.close()
    path = path.replace(".c",".h")
    f = open(path,"w+")
    s = "#ifndef "+os.path.basename(path).replace(".h","_h_")+"\n"
    s += "#define "+os.path.basename(path).replace(".h","_h_")+"\n"
    s += "\n#include <fixpm.h>\n\n"
    for l,x in enumerate(table.keys()):
        if l == len(table.keys())-1:
            s+="sg32 "+x+"["+str(len(table[x]))+"];"
        else:
            s+="sg32 "+x+"["+str(len(table[x]))+"];\n"
    f.write(s)
    f.write("\n\n#endif")
    f.close()
    print "\nFiles created: "+name+".c, "+name+".h"

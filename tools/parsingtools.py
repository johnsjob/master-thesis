from serial.tools import list_ports
#--------------------------------------------------------------------------#
import numpy
from numpy import pi as pi
from numpy import array as mat
#--------------------------------------------------------------------------#
import shlex

__file__ = "parsingtools.py"
# Should use the files from submodule
import sys,  os
filePath, fileName=os.path.split(__file__) #you are here
sys.path.append(os.path.normpath(os.path.join(filePath, '../int/penfirmware-zeus/source/test/test_script/')))
from apranalysis import *
from FrameInfo import *
#==========================================================================#
def parse_file(path):
		f = open(path)
		table = read_file_values(f)
		f.close()
		return table
#----------------------------#
def read_file_values(f = None, dtype=float):
		if f is None:
				raise IOError("file-handle invalid.")
		table = {}
		while True:
				line = f.readline()
				if len(line) == 0:
					break

				result = line.split("=")
				if len(result) == 1 and (('-' in result[0] or '#' in result[0] or '\n' in result[0]) or
										 ('[' in result[0] or ']' in result[0])):
					continue
				elif len(result) > 2:
					continue
				elif len(result) == 2:
					name,weight = result
				else:
					name = result                    
				name = name.strip()
				weight = weight.strip()

				#If it is a tuple value...
				if '(' in weight or ')' in weight:
					weight = weight.replace("(","")
					weight = weight.replace(")","")
					weight = weight.split(",")
					for i,w in enumerate(weight):
						# Some line values contains spaces between separators for some reason
						# these are removed if they exist
						w = w.strip()
						# In case a string contains descriptive text
						# at the end of the values we extract only valid values
						w = w.split(' ')[0]
						weight[i] = dtype(w)
				elif '{' in weight or '}' in weight: # Descriptive text more like a title for
													 # a subsection
						continue
				elif '[' in weight or ']' in weight:
					continue
				else:
						weight = dtype(weight)
				if not table.has_key(name):
					table[name] = []
				#table[name].insert(len(table[name])-1,weight)
				table[name].append(weight)
		return table
#----------------------------#
def read_file_fields(f = None):
		if f is None:
				raise IOError("file-linedata invalid.")
		table = {}
		while True:
				line = f.readline()
				if len(line) == 0:
					break
				result = line.split("=")
				if len(result) != 2:
						continue
				name, weight = result
				name = name.strip()
				if not table.has_key(name):
					table[name] = None
		for i,k in enumerate(table):
				table[k] = k
		return table
#==========================================================================#
def list_to_file(name,valuelist,filename='tmp.txt'):
	try:
		if (type(valuelist) not in [list, numpy.ndarray, tuple]
			or str(name) == 'None' or len(name) < 3 or (type(name) is not str)
			or str(filename) == 'None' or len(filename) < 4 or (type(filename) is not str)):
			raise(TypeError("\tFunction only accepts lists, numpy.ndarrays or tuples."\
							"\n\tName_of_data must be minimum 3 characters long."\
							"\n\tFilename can't be empty or none type and must be\n\tminimum 4 characters excluding filetype extension."\
							"\n\tname_of_data and filename must be strings."))
		if(len(filename.split(".")) != 2): #No split occured or invalid split occured
			raise(TypeError("\tInvalid format of filename."))
	except TypeError as e:
			if str(filename) == 'None':
				filename = '<NO FILENAME GIVEN>'
			print "Usage:\n\tlist_to_file(<name_of_data>, <list_of_values>)"
			print "Description:\n\tWrites raw-data to "+filename
			print "Returns: Number of charatcers (bytes) written."
			print "\nPython Error:\n"+str(e)
	ret = 0;
	f = open(filename,"w+")
	for i,v in enumerate(valuelist):
		if i != len(valuelist)-1:
			data = name+" = "+str(v)+"\n"
		else:
			data = name+" = "+str(v)
		f.write(data)
		ret += len(data)
	f.close()
	return ret
#--------------------------------------------------------------------------#
def load_list_of_matrices(filename):
	#Reads the file
	f = open(filename, "r")
	lines = f.readlines()
	f.close()

	#removes newlines
	for line in lines:
		line = line.strip()

	#removes elements containing '========'
	filtered = [x for x in lines if any(c.isdigit() for c in x)]

	#Each element in filtered is a string of digits separated by spaces
	#This converts these strings to lists of digits
	newlist = list()
	for l in filtered:
		subl = shlex.split(l)
		subl = [float(i) for i in subl]
		#print subl
		newlist.append(subl)

	#Creates a new list of numpy matrices from newlist
	res = []
	for i in range(0,len(newlist),3):
		e = [newlist[i],newlist[i+1],newlist[i+2]]
		if len(e) < 3:
			continue
		try:
			e = numpy.reshape(e,(3,3))
		except:
			print len(e)
			print e
		#print e
		res.append(e)

	#Return the result
	return res
#--------------------------------------------------------------------------#
def get_input_path(filename):
	inputpath = "../resources/pen_calculated_orientation_data/"
	typ = ".txt"
	name  = filename
	return inputpath + name + typ
#--------------------------------------------------------------------------#
def get_fixedpoint_input_path(filename):
	inputpath = "../resources/pen_calculated_orientation_data_fixedpoint/"
	typ = ".txt"
	name  = filename
	return inputpath + name + typ
#--------------------------------------------------------------------------#
def get_rawdata_path(filename):
	inputpath = "../resources/pen_T-matrix_measurements/"
	typ = ".txt"
	name  = filename
	return inputpath + name + typ
#--------------------------------------------------------------------------#
def get_testcase_path(filename):
	inputpath = "../resources/test-cases/"
	typ = ".txt"
	name  = filename
	return inputpath + name + typ
#--------------------------------------------------------------------------#
def get_output_error_path(filename):
	inputpath = "../resources/pen_error_data/"
	typ = ".txt"
	name  = filename
	return inputpath + name + typ
#--------------------------------------------------------------------------#
def get_fixedpoint_output_error_path(filename):
	inputpath = "../resources/pen_error_data_fixedpoint/"
	typ = ".txt"
	name  = filename
	return inputpath + name + typ
#--------------------------------------------------------------------------#
def get_substr_elements(subl):
	return subl.split(" ")
#--------------------------------------------------------------------------#
def remove_invalid_elements(subl):
	return [x for x in subl if any(c.isdigit() for c in x)]
#--------------------------------------------------------------------------#
def remove_invalid_elements2(subl):
	return [x for x in subl if any(c.isdigit() for c in x) and "." in x]
#--------------------------------------------------------------------------#
def get_elements_floats(subl):
	l = [x for x in subl if any(c.isdigit() for c in x)]
	return [float(i) for i in l]
#--------------------------------------------------------------------------#
def get_elements_ints(subl):
	l = [x for x in subl if any(c.isdigit() for c in x)]
	return [int(i) for i in l]
#--------------------------------------------------------------------------#
def process_file(filename):
	f = open(filename, "r")
	lines = f.readlines()
	f.close()
	for index,line in enumerate(lines):
		line = line.replace(":",""); #remove ':'
		line = line.replace("[",""); #remove '['
		line = line.replace("]",""); #remove ']'
		line = line.replace("(",""); #remove '('
		line = line.replace(")",""); #remove '('
		line = line.replace("#",""); #remove '('
		line = line.strip() #remove newline
		lines[index] = line.replace("=",""); #remove '='
	return lines
#--------------------------------------------------------------------------#
def get_matrix(index,lines):
	#print "EXTRACTED MATRIX: "
	#print index
	l1 = get_substr_elements( lines[index+1] )
	l1 = get_elements_floats( l1 )
	#print l
	l2 = get_substr_elements( lines[index+2] )
	l2 = get_elements_floats( l2 )
	#print l
	l3 = get_substr_elements( lines[index+3] )
	l3 = get_elements_floats( l3 )
	#print mat([l1,l2,l3])

	ret = mat([l1,l2,l3])
	return ret
#--------------------------------------------------------------------------#
def get_matrix_int(index,lines):
	#print "EXTRACTED MATRIX: "
	#print index
	l1 = get_substr_elements( lines[index+1] )
	l1 = get_elements_ints( l1 )
	#print l
	l2 = get_substr_elements( lines[index+2] )
	l2 = get_elements_ints( l2 )
	#print l
	l3 = get_substr_elements( lines[index+3] )
	l3 = get_elements_ints( l3 )
	#print mat([l1,l2,l3])

	ret = mat([l1,l2,l3])
	return ret
#--------------------------------------------------------------------------#
def get_values(subl,fieldname):
	for index,l in enumerate(subl):
		if fieldname in l:
			#print "FIELDNAME:"+fieldname
			ret = list(subl[index:len(subl)])
			#print "PRE-RET: "+str(ret)
			tmp = list()
			for x in ret:
				#print "UNSPLIT: "+ str(x)
				y = x.split(",")
				#print "SPLIT: "+ str(y)
				y = remove_invalid_elements2(y)
				#print "SPLIT VALID: "+ str(y)
				#print "SPLIT VALID LEN: "+ str(len(y))
				if len(y) > 0:
					for z in y:
						tmp.append(float(z))
					#print "TMP:" + str(tmp)
					if len(tmp)>1:
						ret = tmp
					else:
						ret = tmp[0]
					#print "DATA: " + str(ret)
			return ret
	return None
#--------------------------------------------------------------------------#
def get_values_int(subl,fieldname):
	for index,l in enumerate(subl):
		if fieldname in l:
			#print "FIELDNAME:"+fieldname
			ret = list(subl[index:len(subl)])
			#print "PRE-RET: "+str(ret)
			tmp = list()
			for x in ret:
				#print "UNSPLIT: "+ str(x)
				y = x.split(",")
				#print "SPLIT: "+ str(y)
				y = remove_invalid_elements(y)
				#print "SPLIT VALID: "+ str(y)
				#print "SPLIT VALID LEN: "+ str(len(y))
				if len(y) > 0:
					for z in y:
						tmp.append(int(z))
					#print "TMP:" + str(tmp)
					if len(tmp)>1:
						ret = tmp
					else:
						ret = tmp[0]
					#print "DATA: " + str(ret)
			return ret
	return None
#--------------------------------------------------------------------------#
def get_float_from_field(lines,fieldname):
	res = []
	for index,line in enumerate(lines):
		subl = get_substr_elements(line)
		#print "SUBL: "+str(subl)
		ret = get_values(subl,fieldname)
		if ret is not None:
			#print "RET: "+str(ret)
			res.append(ret)
	return res
#--------------------------------------------------------------------------#
def get_int_from_field(lines,fieldname):
	res = []
	for index,line in enumerate(lines):
		subl = get_substr_elements(line)
		#print "SUBL: "+str(subl)
		ret = get_values_int(subl,fieldname)
		if ret is not None:
			#print "RET: "+str(ret)
			res.append(ret)
	return res
#--------------------------------------------------------------------------#
def get_matrix_from_field(lines,fieldname):
	res = []
	for index,line in enumerate(lines):
		subl = get_substr_elements(line)
		#print subl
		if( fieldname in subl[0]):
			#print fieldname+"-MATRIX EXTRACTED!"
			#print get_matrix(index,lines)
			ret = get_matrix(index,lines)
			res.append(ret)
	return res
#--------------------------------------------------------------------------#
def get_matrix_from_field_int(lines,fieldname):
	res = []
	for index,line in enumerate(lines):
		subl = get_substr_elements(line)
		#print subl
		if( fieldname in subl[0]):
			#print fieldname+"-MATRIX EXTRACTED!"
			#print get_matrix(index,lines)
			ret = get_matrix_int(index,lines)
			res.append(ret)
	return res
#--------------------------------------------------------------------------#
def get_all_angles(lines):
	res = []
	tmp = []
	for index,line in enumerate(lines):
		subl = get_substr_elements(line)
		#print subl
		#print "EXTRACTED VALUE(s):"
		ret = get_values(subl,"tilt")
		ret2 = get_values(subl,"skew")
		ret3 = get_values(subl,"rot")
		if ret is not None:
#            print "RET: "+str(ret)
			tmp.append(ret)
		if ret2 is not None:
#            print "RET2: "+str(ret2)
			tmp.append(ret2)
		if ret3 is not None:
#            print "RET3: "+str(ret3)
			tmp.append(ret3)
#            print "TMP: "+str(tmp)
			res.append(tmp)
			#print "RET: "+str(tmp)
			#print tmp
			tmp = []
	return res
#--------------------------------------------------------------------------#
def get_all_T(lines):
	l_T = get_matrix_from_field(lines,"T-MAT")
	return l_T
#--------------------------------------------------------------------------#
def get_all_R(lines):
	l_R = get_matrix_from_field(lines,"R-MAT")
	return l_R
#--------------------------------------------------------------------------#
def get_all_R_int(lines):
	l_R = get_matrix_from_field_int(lines,"R-MAT")
	return l_R
#--------------------------------------------------------------------------#
def get_all_Mv(lines):
	l_mag = get_float_from_field(lines,"Mv")
	return l_mag
#--------------------------------------------------------------------------#
def get_all_xyro(lines):
	l_xyro = get_float_from_field(lines,"xro,yro")
	return l_xyro
#--------------------------------------------------------------------------#
def load_list_of_matrices(filename):
	#Reads the file
	f = open(filename, "r")
	lines = f.readlines()
	f.close()

	#removes newlines
	for line in lines:
		line = line.strip()

	#removes elements containing '========'
	filtered = [x for x in lines if any(c.isdigit() for c in x)]

	#Each element in filtered is a string of digits separated by spaces
	#This converts these strings to lists of digits
	newlist = list()
	for l in filtered:
		subl = shlex.split(l)
		subl = [float(i) for i in subl]
		#print subl
		newlist.append(subl)

	#Creates a new list of numpy matrices from newlist
	res = []
	for i in range(0,len(newlist),3):
		e = [newlist[i],newlist[i+1],newlist[i+2]]
		if len(e) < 3:
			continue
		try:
			e = numpy.reshape(e,(3,3))
		except:
			print len(e)
			print e
		#print e
		res.append(e)

	#Return the result
	#print str(type(res))
	return res
#--------------------------------------------------------------------------#
def pr(st,ind=0):
        s = ""
        for k in range(ind):
                s += "\t"
        out = s+st
        c = 80
        if len(out) > c:
                print out[0:c]+" ... "
                pr(" ... "+out[c:], ind)
        else:
                print out
#--------------------------------------------------------------------------#
def newl():
        print ""
#--------------------------------------------------------------------------#
if __name__ == "__main__":
	pr("=== Running legend code ===")
	newl()
	
	filepath = "../resources/pen_calculated_orientation_data/case0__theta_is_largerThan_45__skew_is_inRange_0_to_minus_90__rot_is_near_0.txt"
	pr("file_path = " + filepath,1)
	
	lines = process_file(filepath)
	list_of_T = get_all_T(lines)

	newl()
	pr("loaded matrices =  " + str(len(list_of_T)),1)

	newl()
	pr("=== Running newer code ===")
	table = parse_file(filepath)

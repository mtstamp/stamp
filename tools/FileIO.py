#############################################################################
#
#Updated by: Yiqin Wang
#
##############################################################################

import sys
import subprocess
import types
import struct
import gzip,bz2,zipfile

def execute(cmd, stderr = sys.stderr, run = True):
	#a function wrapper to run subprocess.call
	print >>stderr, "CMD: %s" % cmd
	if (run):
		return subprocess.call(cmd, shell=True)
	else:
		return None

def pipe_input(cmd, stderr = sys.stderr, run = True):
	#a function wrapper to run subprocess.Popen for stdin
	print >>stderr, "CMD: %s" % cmd
	if (run):
		pf = subprocess.Popen(cmd, stdin=subprocess.PIPE, stderr=None, shell=True)
		return pf
	else:
		return None
	
def pipe_output(cmd, stderr = sys.stderr, run = True):
	#a function wrapper to run subprocess.Popen for stdout
	print >>stderr, "CMD: %s" % cmd
	if (run):
		pf = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=None, shell=True)
		return pf
	else:
		return None

def read_arguments(arg_file):
	arguments = {}
	fp = FileIO(arg_file, "r")
	for line in fp:
		line = line.rstrip("\r\n")
		if (not line):
			continue
		if (line[0] == "#"):
			continue
		n = line.find("=")
		assert n >=0, "Unknown argument line: %s" % line
		arg = line[:n].strip()
		value = line[(n+1):].strip()
		arguments[arg] = value
	fp.close()
	return arguments

def write_arguments(arguments, arg_file):
	fp = FileIO(arg_file, "w")
	for i in arguments:
		if (len(i) == 1):
			arg = i[0]
			value = ""
			comment = ""
		elif (len(i) == 2):
			arg = i[0]
			value = i[1]
			comment = ""
		elif (len(i) >= 3):
			arg, value, comment = i[:3]
		if (comment):
			fp.write("#" + comment + "\n")
		if (arg):
			if (value is None):
				value = ""
			fp.write(arg + "=" + str(value) + "\n")
	fp.close()

class FileIO:
	""" 
	This is a class for handling various types of file input and output
	
	Attributes
	----------
	f: file handle (text or compressed files gz/bz/zip)
	random_access: whether supports random access in the file
	
	"""
	def __init__(self, fname, mode, type = None, **kargs):
		if (type == "gz" or fname.endswith(".gz")):
			self.f = gzip.open(fname, mode, **kargs)
			self.random_access = True
		elif (type == "bz" or fname.endswith(".bz2")):
			self.f = bz2.BZ2File(fname, mode, **kargs)
			self.random_access = True
		elif (type == "zip" or fname.endswith(".zip")):
			self.f = zipfile.ZipFile(fname, mode, **kargs)
			self.random_access = False
		else:
			self.f = open(fname, mode)
			self.random_access = True
		self.name = fname
	
	def write(self, buffer):
		#a function wrapper for write
		return self.f.write(buffer)
	
	def read(self, length):
		#a function wrapper for read
		return self.f.read(length)
	
	def seek(self, offset, whence=0):
		#a function wrapper for seek
		if (self.random_access):
			return self.f.seek(offset, whence) #not support zip file
		else:
			return -1
	
	def tell(self):
		#a function wrapper for tell
		if (self.random_access):
			return self.f.tell() #not support zip file
		else:
			return -1
	
	def close(self):
		#a function wrapper for close
		return self.f.close()
	
	def readline(self):
		#a function wrapper for readline
		return self.f.readline()
		
	def __iter__(self):
		#a function wrapper for __iter__
		return self
	
	def next(self):
		#a function wrapper for next
		line = self.readline()
		if (len(line) == 0):
				raise StopIteration
		return line

class BinIO(FileIO):
	""" 
	This is a class for writing and reading binary data in a file
	
	Attributes
	----------
	endian: big, small or system-based 
	other inherited attributes from FileIO
	
	"""
	def __init__(self, fname, mode, endian):
		assert endian in ["=", "<", ">"], "Unknown endian"
		self.__endian__ = endian
		FileIO.__init__(self, fname, mode)
	
	def getDigits(self, ctype, line, offset = 0, num = 1):
		'''
		read digits from line
		'''
		if (ctype in ['c','b','B','?']):
			count = 1
		elif (ctype in ['h','H']):
			count = 2
		elif (ctype in ['i','I','l','L','f']):
			count = 4
		elif (ctype in ['q','Q','d']):
			count = 8
		else:
			return None
		fmt = self.__endian__ + ctype*num
		x = struct.unpack_from(fmt,line,offset=offset)
		if (num == 1):
			return x[0]
		else:
			return x
	
	def readDigits(self, ctype, num = 1):
		'''
		read digits from current place in file
		'''
		if (ctype in ['c','b','B','?']):
			count = 1
		elif (ctype in ['h','H']):
			count = 2
		elif (ctype in ['i','I','l','L','f']):
			count = 4
		elif (ctype in ['q','Q','d']):
			count = 8
		else:
			return None
		fmt = self.__endian__ + ctype*num
		line = self.read(count*num)
		x = struct.unpack_from(fmt,line)
		if (num == 1):
			return x[0]
		else:
			return x
	
	def putDigits(self, ctype, x, line = None, offset = 0, num = 1):
		'''
		write digits to line
		'''
		if (ctype in ['c','b','B','?']):
			count = 1
		elif (ctype in ['h','H']):
			count = 2
		elif (ctype in ['i','I','l','L','f']):
			count = 4
		elif (ctype in ['q','Q','d']):
			count = 8
		else:
			return None
		fmt = self.__endian__ + ctype*num
		if (not isinstance(x, types.ListType) and not isinstance(x, types.TupleType)):
			x = [x,]
		if (line is None):
			line = struct.pack(fmt, *x)
		else:
			struct.pack_into(fmt, line, offset, *x)
		return line
	
	def writeDigits(self, ctype, x, num):
		'''
		write digits to current place in file
		'''
		if (ctype in ['c','b','B','?']):
			count = 1
		elif (ctype in ['h','H']):
			count = 2
		elif (ctype in ['i','I','l','L','f']):
			count = 4
		elif (ctype in ['q','Q','d']):
			count = 8
		else:
			return None
		fmt = self.__endian__ + ctype*num
		if (not isinstance(x, types.ListType) and not isinstance(x, types.TupleType)):
			x = [x,]
		line = struct.pack(fmt, *x)
		self.write(line)
		return line
	
	def getByte(self, line, offset=0, num=1):
		return self.getDigits('b', line, offset, num)
	def getUByte(self, line, offset=0, num=1):
		return self.getDigits('B', line, offset, num)
	def getShort(self, line, offset=0, num=1):
		return self.getDigits('h', line, offset, num)
	def getUShort(self, line, offset=0, num=1):
		return self.getDigits('H', line, offset, num)
	def getInt(self, line, offset=0, num=1):
		return self.getDigits('i', line, offset, num)
	def getUInt(self, line, offset=0, num=1):
		return self.getDigits('I', line, offset, num)
	def getInt64(self, line, offset=0, num=1):
		return self.getDigits('q', line, offset, num)
	def getUInt64(self, line, offset=0, num=1):
		return self.getDigits('Q', line, offset, num)
	def getFloat(self, line, offset=0, num=1):
		return self.getDigits('f', line, offset, num)
	def getDouble(self, line, offset=0, num=1):
		return self.getDigits('d', line, offset, num)
	
	def readByte(self, num=1):
		return self.readDigits('b', num)
	def readUByte(self, num=1):
		return self.readDigits('B', num)
	def readShort(self, num=1):
		return self.readDigits('h', num)
	def readUShort(self, num=1):
		return self.readDigits('H', num)
	def readInt(self, num=1):
		return self.readDigits('i', num)
	def readUInt(self, num=1):
		return self.readDigits('I', num)
	def readInt64(self, num=1):
		return self.readDigits('q', num)
	def readUInt64(self, num=1):
		return self.readDigits('Q', num)
	def readFloat(self, num=1):
		return self.readDigits('f', num)
	def readDouble(self, num=1):
		return self.readDigits('d', num)
	
	def putByte(self, x, line = None, offset = 0, num=1):
		return self.putDigits('b', x, line, offset, num)
	def putUByte(self, x, line = None, offset = 0, num=1):
		return self.putDigits('B', x, line, offset, num)
	def putShort(self, x, line = None, offset = 0, num=1):
		return self.putDigits('h', x, line, offset, num)
	def putUShort(self, x, line = None, offset = 0, num=1):
		return self.putDigits('H', x, line, offset, num)
	def putInt(self, x, line = None, offset = 0, num=1):
		return self.putDigits('i', x, line, offset, num)
	def putUInt(self, x, line = None, offset = 0, num=1):
		return self.putDigits('I', x, line, offset, num)
	def putInt64(self, x, line = None, offset = 0, num=1):
		return self.putDigits('q', x, line, offset, num)
	def putUInt64(self, x, line = None, offset = 0, num=1):
		return self.putDigits('Q', x, line, offset, num)
	def putFloat(self, x, line = None, offset = 0, num=1):
		return self.putDigits('f', x, line, offset, num)
	def putDouble(self, x, line = None, offset = 0, num=1):
		return self.putDigits('d', x, line, offset, num)
	
	def writeByte(self, x, num = 1):
		return self.writeDigits('b', x, num)
	def writeUByte(self, x, num = 1):
		return self.writeDigits('B', x, num)
	def writeShort(self, x, num = 1):
		return self.writeDigits('h', x, num)
	def writeUShort(self, x, num = 1):
		return self.writeDigits('H', x, num)
	def writeInt(self, x, num = 1):
		return self.writeDigits('i', x, num)
	def writeUInt(self, x, num = 1):
		return self.writeDigits('I', x, num)
	def writeInt64(self, x, num = 1):
		return self.writeDigits('q', x, num)
	def writeUInt64(self, x, num = 1):
		return self.writeDigits('Q', x, num)
	def writeFloat(self, x, num = 1):
		return self.writeDigits('f', x, num)
	def writeDouble(self, x, num = 1):
		return self.writeDigits('d', x, num)
	
	def readString(self, count):
		line = self.read(count)
		t = count * "s"
		x = "".join(list(struct.unpack(t, line)))
		return x
	
	def readLongString(self):
		string_len = self.readUInt()
		return self.readString(string_len)
	
	def readShortString(self):
		string_len = self.readUShort()
		return self.readString(string_len)
	
	def readByteString(self):
		string_len = self.readUByte()
		return self.readString(string_len)
	
	def writeString(self, x):
		count = len(x)
		t = count * "s"
		line = struct.pack(t, x)
		self.write(line)
	
	def writeLongString(self, x):
		string_len = len(x)
		self.writeUInt(x)
		return self.writeString(x)
	
	def writeShortString(self, x):
		string_len = len(x)
		self.writeUShort(x)
		return self.writeString(x)
	
	def writeByteString(self, x):
		string_len = len(x)
		self.writeUByte(x)
		return self.writeString(x)


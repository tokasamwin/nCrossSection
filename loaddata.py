load=False
try:
	from openmc.data import IncidentNeutron as NEUTRON
	mc=True
except:
	mc=False
	print('openMC not installed, falling back to pyENDF6')
import pyENDF6.ENDF6 as endf

import os, sys, numpy as np, matplotlib.pyplot as plt
from scipy.interpolate import interp1d as interp
for p in sys.path:
	potential=os.path.join(p,'nCrossSection')
	if os.path.exists(potential):
		ncspath=potential
		break

def set_path(path):
	'''
	Sets a path variable within the current installation, should be per-user
	'''
	with open(os.path.join(ncspath,'user.path'),'w+') as f:
		f.write(path)

def read_path():
	with open(os.path.join(ncspath,'user.path')) as f:
		neutronicspath=f.read().strip()
	return neutronicspath

class nuclear_directory(object):
	def __init__(self,path):
		'''
		Data object to store file index for nuclear data
		Default data types are endf
		'''
		self.path=path
		self.list=os.listdir(path)
		count=len(list(filter(lambda i: self.list[0] in i,self.list)))
		if '.' in self.list[0]:
			self.format=self.list[0].split('.')[-1] #captures hdf5 and endf
		else:
			self.format='ace'
		self.describe()
		self.get_data=lambda ID:loaddata.get(ID)
	
	def describe(self):
		self.index={}
		for i in self.list:
			getattr(self,self.format)(i)
	
	def endf(self,i):
		var=i.split('.')[0].replace('n-','')
		Z,s,A=var.split('_')
		try:
			Z,A=int(Z),int(A)
			A=int(A)
			m=0
			self.index[(Z,A)]=os.path.join(self.path,i)
			self.index['{}{}'.format(s,A)]=self.index[(Z,A)]
		except:
			Z=int(Z)
			A,m=[int(i) for i in A.split('m')]
			self.index[(Z,A,m)]=os.path.join(self.path,i)
			self.index['{}{}_{}'.format(s,A,m)]=self.index[(Z,A,m)]

	def ace(self,i,post=''):
		if '.' in i: return
		Z=int(i[:2])
		s=i[2:4].replace('_','')
		A=int(i[-3:])
		self.index[(Z,A)]=os.path.join(self.path,i+post)
		self.index['{}{}'.format(s,A)]=self.index[(Z,A)]
		
	def hdf5(self,i):
		name=i.split('.')[0]
		self.ace(name,post='.hdf5')


try:
	neutronicspath=read_path()
	data=nuclear_directory(neutronicspath)
except:
	rtext='Warning: Data unavailable\nSet a path to a neutronics data file'
	print(rtext)
	print(neutronicspath)
	print(nuclear_directory(neutronicspath))


class mcdata(object):
	def __init__(self,ID):
		self.ID=ID
		self.data=getattr(NEUTRON,'from_'+data.format)(data.index[ID])
		self.XStype=self.data.reactions.keys()
	
	def read(self,E,MT=1):
		T=list(self.data[MT].xs.keys())[0]
		return self.data[MT].xs[T](E)

class endfdata(object):
	def __init__(self,ID):
		with open(data.index[ID]) as f:
			self.lines=f.readlines()
		content=endf.list_content(self.lines)
		self.XStype=[i[-1] for i in filter(lambda i: i[1]==3,content)]
		self.gen_arrays()

	def gen_arrays(self):
		self.rawdata={}
		self.data={}
		Edebug=10**np.linspace(-3,0)
		f=plt.figure()
		ax=f.subplots(2)
		for xs in self.XStype:
			self.rawdata[xs]=endf.read_table(endf.find_section(self.lines,MF=3,MT=xs))
			self.data[xs] = interp(*self.rawdata[xs],kind='linear',bounds_error=0,fill_value=0) # construct linear interpreter functions for each cross section
			if self.debug:
				ax[0].plot(*self.rawdata[xs])
				ax[0].set_xlim([0,1])
				ax[0].set_ylim([0,2.5])
				ax[1].plot(Edebug,self.data[xs](Edebug))
				ax[1].set_xlim([0,1])
				ax[1].set_ylim([0,2.5])
			#this step currently doesn't work, don't know how to make it work
			#seems to be an error in how the lambda functions are set up
	
	def read(self,E,MT=1):
		return self.data[MT](E)

if mc:
	def get_data(ID):
		return mcdata(ID)
else:
	def get_data(ID):
		return endfdata(ID)
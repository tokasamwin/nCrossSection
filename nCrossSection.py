#from nova import ENDF6
import numpy as np, os, sys, matplotlib.pyplot as plt, mendeleev
from math import isnan, isinf
amu=1.660539040e-27 # atomic mass unit, used to find number densities
load=True
if load:
	try:
		from openmc.data import IncidentNeutron as NEUTRON
	except:
		print('Cannot load OpenMC, cannot read cross section data')

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
		self.get_data=lambda ID:self.index[ID]
	
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

if load:
	try:
		neutronicspath=read_path()
		data=nuclear_directory(neutronicspath)
	except:
		rtext='Warning: Data unavailable\nSet a path to a neutronics data file'
		print(rtext)
		print(neutronicspath)
		print(nuclear_directory(neutronicspath))

class isotope(object):
	def __init__(self,Z,A,m=None):
		self.Z=Z
		self.A=A
		self.mnum=m
		self.id=join(self.Z,self.A,self.mnum)
		self.m=A*amu
		self.notloaded=True
		if load:
			self.read()
			
	def read(self):
		#finding cross sections available
		self.data=getattr(NEUTRON,'from_'+data.format)(data.index[self.id])
		self.XStype=self.data.reactions.keys()
		self.notloaded=False
	
	def XSplot(self,*MT,**kwargss):
		pyplot_objs=[]
		if 'fig' in kwargs:
			f=kwargs.get('fig')
			pyplot_objs.append(f)
			if len(f.axes):
				ax=f.axes[0]
			else:
				ax=f.subplots()
			pyplot_objs.append(ax)
		elif 'ax' in kwargs:
			ax=kwargs.get('ax')
			pyplot_objs.append(ax)
		else:
			f=plt.figure()
			ax=f.subplots()
			pyplot_objs+=[f,ax]
		if self.notloaded:
			self.read()
		for i in MT:
			try:
				ax.loglog(self.E[i],self.XSdata[i],label='{},{}:{}'.format(self.Z,self.A,i))
			except:
				...
		if not len(ax.get_xlabel()):
			ax.xlabel('Energy ($eV$)')
			ax.ylabel('Microscopic cross section (barn or $b$)')
			ax.legend(loc='best')
		for kw in kwargs:
			for obj in pyplot_objs:
				if hasattr(obj,kw):
					getattr(obj,kw)(kwargs[kw])
				elif hasattr(obj,'set_'+kw):
					getattr(obj,kw)(kwargs[kw])
	
	def XSfind(self,E,MT=1):
		if self.notloaded:
			self.read()
		T=list(self.data[MT].xs.keys())[0]
		return self.data[MT].xs[T](E)
	
class element:
	def __init__(self,Z,*A_data):
		'''
		Always returns microscopic cross sections of mixtures
		Multiply by bulk atomic density to determine macroscopic cross section of mixtures
		
		Inputs:
		Z is atomic number (proton number)
		A_tuples is a tuple of (A,c) where c is the composition of A
			A is isotopic mass
			c is the composition of isotopic mass
			composition is normalised, so can specify as (6,0.25),(7,0.75) or (6,1),(7,3)
		if tuples aren't given, the element defaults to typical composition (see 'mendeleev' package)
		'''
		if type(Z)==int:			
			self.Z=Z
		elif type(Z)==str:
			self.Z=mendeleev.element(Z).atomic_number
			Z=self.Z
		if len(A_data): #equivalent to len(A_data)>0, but faster
			A,comp=zip(*A_data)
		else:
			'''
			This is the defaulting behaviour with no data provided:
				revert to standard isotope composition
			'''
			mend=mendeleev.element(Z)
			zipped=[(iso.mass_number,iso.abundance) for iso in mend.isotopes]
			filtered=filter(lambda i:i[-1]!=None,zipped)
			A,comp=zip(*filtered)
		comp=[i/sum(comp) for i in comp]
		self.dataav=[]
		self.avm=0
		self.iso={}
		self.comp={}
		for i,x in zip(A,comp):
			if (Z,i) in data.index:
				self.iso[(Z,i)]=isotope(Z,i)
				self.comp[(Z,i)]=x
				self.avm+=self.iso[(Z,i)].m*self.comp[(Z,i)]
				if load:
					for MT in self.iso[(Z,i)].XStype:
						if MT not in self.dataav:
							self.dataav.append(MT)
		self.dataav.sort()

class compound:
	def __init__(self,rho,*e_data,mix=None):
		importconds=all([len(e_data)==0,rho==None,mix is not None])
		if not importconds:
			try:
				self.defaultinit(rho,*e_data)
			except:
				raise AttributeError('Need input of elements, composition array and density of compound, or a mixture!')
		elif importconds:
			self.importfrommix(mix)
		else:
			print('Something failed')
	
	def defaultinit(self,rho,*e_data):
		'''
		Inputs:
		rho is density in kg/m3
		elemarray is an element object or list of element objects
			These element objects must be initialised before generating a compound
		elemcomp is the composition fraction of each item in elemarray,
			which can be either a percentage of total, fraction 
			or as in a chemical formula (i.e. H2O would be 2,1, U3O8 would be 3,8)
			the code takes the elemarray and normalises it for appropriate ratios
		e_data is a set of tuples
		'''
		elems,comp=zip(*e_data)
		comp=[i/sum(comp) for i in comp]
		self.rho=rho
		self.species={}
		self.dataav=[]
		for e,c in zip(elems,comp):
			for i in e.iso:
				self.species[i]=e.iso[i]
			for MT in list(e.dataav):
				if MT not in self.dataav:
					self.dataav.append(MT)
		self.dataav.sort() 
		#finding the isotope composition in a compound
		self.isocomp={}
		for e,c in zip(elems,comp):
			for iso in e.iso:
				if iso in self.isocomp:
					self.isocomp[iso]+=c*e.comp[iso]
				else:
					self.isocomp[iso]=c*e.comp[iso]
		totalcomp=sum(self.isocomp.values())
		self.isofrac={}
		for iso in self.isocomp:
			self.isofrac[iso]=self.isocomp[iso]/totalcomp
		fracsum=0
		for (Z,A) in self.isocomp:
			fracsum+=A*self.isofrac[(Z,A)]
		self.N=self.rho/(amu*fracsum)
		for iso in self.isocomp:
			self.isocomp[iso]=self.isofrac[iso]*self.N
	
	def totalXS(self,E,MT=1):
		totXS=0
		for x,iso in zip(self.isofrac,self.species):
			try:
				totXS+=self.N*self.isofrac[x]*self.species[iso].XSfind(E,MT=MT)/10**28
			except:
				totXS+=0
		return totXS
	
	def importfrommix(self,mix):
		#import number density
		self.N=mix.N
		#import data available
		self.dataav=mix.dataav
		self.species=mix.species
		self.rho=mix.rho
		self.isocomp=mix.isocomp
		self.isofrac=mix.isofrac

class mixture:
	def __init__(self,*c_data,Eres=1000,label='Default'):
		'''
		c_data are tuples or iterables of length: 2, holding:
			1: Compound (object, initialised by ncs.compound)
			2: Volume representation (float or int)
		The volume representations are normalised, e.g. norm([0.1,0.9]) = norm([10,90]) = [0.1,0.9]
		For a pebble bed, you may describe the compounds as:
			Pebbles (Be)
			Fill gas (e.g. He)
			Structural support (Eurofer)
		Each of these are compounds, made up of elements
		The mixture is then generated as such:
		pebblebed=nCrossSection.mixture([pebbles,gas,SS],[75,15,10])
		'''
		self.gen=False
		self.label=label
		self.N=0
		
		comps,vols=zip(*c_data)
		self.vol_frac_array=[i/sum(vols) for i in vols]
		self.comp_array=comps
		self.res=Eres
		self.dataav=[]
		self.rho=0
		self.isocomp={}
		self.species={}
		for comp,frac in zip(self.comp_array,self.vol_frac_array):
			self.N+=comp.N*frac
			for MT in comp.dataav:
				if MT not in self.dataav:
					self.dataav.append(MT)
			self.rho+=comp.rho*frac
			for iso in comp.isocomp:
				if iso in self.isocomp:
					self.isocomp[iso]+=frac*comp.isocomp[iso]
				else:
					self.isocomp[iso]=frac*comp.isocomp[iso]
					self.species[iso]=comp.species[iso]
		self.dataav+=[101] if load else []
		self.dataav.sort() if load else ...
		totalN=sum(self.isocomp.values())
		self.isofrac={}
		self.isomass={}
		for iso in self.isocomp:
			self.isofrac[iso]=self.isocomp[iso]/totalN
			self.isomass[iso]=self.isocomp[iso]*amu*iso[1]
		self.XS={}
			
	def mixXS(self,E,MT=1):
		avXS=0
		if MT!=101:
			for iso in self.species:
				if MT in self.species[iso].data:
					avXS+=self.isocomp[iso]*self.species[iso].XSfind(E,MT=MT)/10**28
		else:
			for iso in self.species:
				for MT in range(102,118):
					if MT in self.species[iso].data:
						avXS+=self.isocomp[iso]*self.species[iso].XSfind(E,MT=MT)/10**28
		return avXS
	
	def XSgen(self,*MT):
		lowlim=np.log10(0.01)
		uplim=np.log10(15e6)
		exparr=np.linspace(lowlim,uplim,int(self.res))
		self.Earr=10**exparr
		for mt in MT:
			if mt not in self.dataav:
				continue
			self.XS[mt]=[self.mixXS(e,MT=mt) for e in self.Earr]
				
	def XSplot(self,*MT,disp=['MT','label'],**kwargs):
		pyplot_objs=[]
		if 'fig' in kwargs:
			f=kwargs.get('fig')
			pyplot_objs.append(f)
			if len(f.axes):
				ax=f.axes[0]
			else:
				ax=f.subplots()
			pyplot_objs.append(ax)
		elif 'ax' in kwargs:
			ax=kwargs.get('ax')
			pyplot_objs.append(ax)
		else:
			f=plt.figure()
			ax=f.subplots()
			pyplot_objs+=[f,ax]
		label=''
		for r in MT:
			if 'MT' in disp:
				label+='MT{}; '.format(r)
			if 'label' in disp:
				label+=self.label
			if r not in self.XS:
				self.XSgen(r)
			ax.loglog(self.Earr,self.XS[r],label=label.format(r),lw=0.5)
		if not len(ax.get_xlabel()):
			ax.set_xlabel('Neutron energy ($eV$)')
			ax.set_ylabel('Macroscopic Cross Section ($1/m$)')
			ax.set_xlim([min(self.Earr),max(self.Earr)])
		ax.legend(loc='best')
		for kw in kwargs:
			for obj in pyplot_objs:
				if hasattr(obj,kw):
					getattr(obj,kw)(kwargs[kw])
				elif hasattr(obj,'set_'+kw):
					getattr(obj,kw)(kwargs[kw])

def join(*args):
	'''
	Generic tool for joining non-empty entries into strings
	'''
	args=tuple(filter(lambda i:i!=None,args))
	if type(args[0]) is str:
		return ''.join([str(arg) for arg in args])
	else:
		return args
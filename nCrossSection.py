#from nova import ENDF6
import numpy as np, os, sys, matplotlib.pyplot as plt, mendeleev
from math import isnan, isinf

import nCrossSection.loaddata as ld

amu=1.660538782-27 # atomic mass unit, used to find number densities
load=False
Edef=10**np.linspace(-2,7,num=200)
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
		if self.notloaded: self.data=ld.get_data(self.id)
		self.XSfind=self.data.read
		self.XStype=self.data.XStype
		#self.data=getattr(NEUTRON,'from_'+data.format)(data.index[self.id])
		#self.XStype=self.data.reactions.keys()
		self.notloaded=False
	
	def XSplot(self,*MT,**kwargs):
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
		if not MT:
			MT=[1]
		if self.notloaded:
			self.read()
		for i in MT:
			ax.loglog(Edef,self.XSfind(Edef,MT=i),label='{},{}:{}'.format(self.Z,self.A,i))
		if not len(ax.get_xlabel()):
			ax.set_xlabel('Energy ($eV$)')
			ax.set_ylabel('Microscopic cross section (barn or $b$)')
			ax.legend(loc='best')
		for kw in kwargs:
			for obj in pyplot_objs:
				if hasattr(obj,kw):
					getattr(obj,kw)(kwargs[kw])
				elif hasattr(obj,'set_'+kw):
					getattr(obj,kw)(kwargs[kw])
	
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
			self.iso[(Z,i)]=isotope(Z,i)
			self.comp[(Z,i)]=x
			self.avm+=self.iso[(Z,i)].m*self.comp[(Z,i)]
			if load:
				if (Z,i) in ld.data.index:
					for MT in self.iso[(Z,i)].XStype:
						if MT not in self.dataav:
							self.dataav.append(MT)
		self.dataav.sort()
	
	def XSfind(self,E,MT=1):
		return sum([self.iso[i].XSfind(E,MT=MT) * self.comp[i] for i in self.iso])

class compound:
	def __init__(self,rho,*e_data,mix=None,label=None):
		self.label=label if label else 'default'
		importconds=all([len(e_data)==0,rho==None,mix is not None])
		if not importconds:
			try:
				self.defaultinit(rho,*e_data)
			except:
				raise AttributeError('Need input of density and elements, or a mixture!')
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
	
	def mixXS(self,E,MT=1):
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
				if MT in self.species[iso].data.__iter__:
					avXS+=self.isocomp[iso]*self.species[iso].XSfind(E,MT=MT)/10**28
		else:
			for iso in self.species:
				for MT in range(102,118):
					if MT in self.species[iso].data.__iter__:
						avXS+=self.isocomp[iso]*self.species[iso].XSfind(E,MT=MT)/10**28
		return avXS
	
	def XSgen(self,*MT):
		for mt in MT:
			if mt not in self.dataav:
				continue
			self.XS[mt]=self.mixXS#[self.mixXS(e,MT=mt) for e in self.Earr]
		
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
		if not MT:
			MT=[1]
		for r in MT:
			if 'MT' in disp:
				label+='MT{}; '.format(r)
			if 'label' in disp:
				label+=self.label
			if r not in self.XS:
				self.XSgen(r)
			ax.loglog(Edef,self.XS[r](Edef),label=label.format(r),lw=0.5)
		if not len(ax.get_xlabel()):
			ax.set_xlabel('Neutron energy ($eV$)')
			ax.set_ylabel('Macroscopic Cross Section ($1/m$)')
			ax.set_xlim([min(Edef),max(Edef)])
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


#mend=mendeleev.element(Z)
#zipped=[(iso.mass_number,iso.abundance) for iso in mend.isotopes]
#filtered=filter(lambda i:i[-1],zipped)
#A,comp=zip(*filtered)
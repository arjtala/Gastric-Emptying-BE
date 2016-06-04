import numpy as np
from scipy.integrate import odeint
from sympy.functions.special.delta_functions import Heaviside
import time

# Define custom Heaviside function
def u(t):
	t = np.asarray(t)
	is_scalar = False if t.ndim > 0 else True
	t.shape = (1,)*(1-t.ndim) + t.shape
	unit_step = np.arange(t.shape[0])
	lcv = np.arange(t.shape[0])
	for place in lcv:
		if t[place] == 0:
			unit_step[place] = .5
		elif t[place] > 0:
			unit_step[place] = 1
		elif t[place] < 0:
			unit_step[place] = 0
	return (unit_step if not is_scalar else (unit_step[0] if unit_step else Heaviside(t)))

# System of equations
class massBalance:

	def __init__(self, t, to, Kpel, PermRates, KgeParams, LagParams, vol, M0, title):
		self.title = title

		self.vol = vol

		self.t = t
		self.to = to

		self.M0 = M0

		# Gastric emptying parameters
		self.p150,self.p1200 = KgeParams[0],KgeParams[6]
		self.p250, self.p2200 = KgeParams[1],KgeParams[7]
		self.p350, self.p3200 = KgeParams[2],KgeParams[8]
		self.theta50, self.theta200 = KgeParams[3],KgeParams[9]
		self.tau50, self.tau200 = KgeParams[4],KgeParams[10]
		self.s50, self.s200 = KgeParams[5],KgeParams[11]
		
		# Delay parameters
		self.a50,self.a200 = LagParams[0],LagParams[6]
		self.b50,self.b200 = LagParams[1],LagParams[7]
		self.c50,self.c200 = LagParams[2],LagParams[8]
		self.d50,self.d200 = LagParams[3],LagParams[9]
		self.e50,self.e200 = LagParams[4],LagParams[10]
		self.f50,self.f200 = LagParams[5],LagParams[11]

		self.PermRates = PermRates
		self.perm0,self.perm1,self.perm2,self.perm3,self.perm4,self.perm5,self.perm6 = self.PermRates

		self.KIntFactors = np.array([0.77, 0.76, 0.75, 0.74, 0.73, 0.72, 0.7],float)
		
		self.Kpel = Kpel
			
	# Gastric emptying functions for 50mL and 200mL volumes
	def kge50ml(self,t):
		sum = 0.0
		for k in range(1,26):
			sum += ((-1.0)**k*np.sin(-self.theta50*np.pi*k*(t-self.tau50)))/k+self.p150
		return self.p250*sum**self.p350+self.s50	
	def kge200ml(self,t):
		sum = 0.0
		for k in range(1,26):
			sum += ((-1.0)**k*np.sin(-self.theta200*np.pi*k*(t-self.tau200)))/k+self.p1200
		return self.p2200*sum**self.p3200+self.s200

	# Delay functions for 50mL and 200mL volumes
	def tlag50ml(self,t):
		return self.a50 - self.b50/(self.c50+self.d50*np.exp(-self.e50*np.mod(t,2.0/self.theta50)+self.f50))
	def tlag200ml(self,t):
		return self.a200 - self.b200/(self.c200+self.d200*np.exp(-self.e200*np.mod(t,2.0/self.theta200)+self.f200))			
			
	# Delayed gastric emptying function for small particules and solutions
	def Kge(self,t,to):
		T = np.mod(t+to,120)
		tlag50 = self.tlag50ml(t+to)
		kge50 = self.kge50ml(t+to)
		tlag200 = self.tlag200ml(t+to)
		kge200 = self.kge200ml(t+to)
		U200 = u(T-tlag200)
		U50 = u(T-tlag50)
		if (self.vol>=200): return U200*kge200
		else: return U50*kge50

	# Intestinal transit functions
	def KInt(self,t,to,index):
		tlag50 = self.tlag50ml(t+to)
		kge50 = self.kge50ml(t+to)
		tlag200 = self.tlag200ml(t+to)
		kge200 = self.kge200ml(t+to)
		T = np.mod(t+to,120)
		U = u(T-tlag50)
		return self.KIntFactors[index]*(U*kge50)

	# Backflow functions	
	def Q1(self,t,to):
		KInt1 = self.KInt(t,to,1)
		return .15*KInt1
	def Q2(self,t,to):
		KInt3 = self.KInt(t,to,3)
		return .15*KInt3
	def Q3(self,t,to):
		KInt5 = self.KInt(t,to,5)
		return .15*KInt5
		
	def dF(self,F,t):
		self.F = F
		self.t = t		
		DSolns, DSolnInt0, DSolnInt1, DSolnInt2, DSolnInt3, DSolnInt4, DSolnInt5, DSolnInt6, DSolnPlasma = F[0:9]
		
		Kge = self.Kge(t,self.to)
		KInt0 = self.KInt(t,self.to,0)
		KInt1 = self.KInt(t,self.to,1)
		KInt2 = self.KInt(t,self.to,2)
		KInt3 = self.KInt(t,self.to,3)
		KInt4 = self.KInt(t,self.to,4)
		KInt5 = self.KInt(t,self.to,5)
		Kie = self.KInt(t,self.to,6)
		Q1 = self.Q1(t,self.to)
		Q2 = self.Q2(t,self.to)
		Q3 = self.Q3(t,self.to)
		
		Mnew = np.array([-DSolns*Kge, \
			DSolns*Kge - DSolnInt0*(KInt0 + self.perm0), \
			DSolnInt0*KInt0 - DSolnInt1*(KInt1 + self.perm1) + DSolnInt2*Q1, \
			DSolnInt1*KInt1 - DSolnInt2*(KInt2 + Q1 + self.perm2), \
			DSolnInt2*KInt2 - DSolnInt3*(KInt3 + self.perm3) + DSolnInt4*Q2, \
			DSolnInt3*KInt3 - DSolnInt4*(KInt4 + Q2 + self.perm4), \
			DSolnInt4*KInt4 - DSolnInt5*(KInt5 + self.perm5) + DSolnInt6*Q3, \
			DSolnInt5*KInt5 - DSolnInt6*(Kie + Q3 + self.perm6), \
			DSolnInt0*self.perm0 + DSolnInt1*self.perm1 + DSolnInt2*self.perm2 + DSolnInt3*self.perm3 + DSolnInt4*self.perm4 + DSolnInt5*self.perm5 + DSolnInt6*self.perm6 - DSolnPlasma*self.Kpel],dtype='float')
		return Mnew

	def Fjac(self,F,t):
		self.t = t
		Fmat = np.array([[Kge,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], \
			[Kge, -(KInt0 + self.perm0),0.0,0.0,0.0,0.0,0.0,0.0,0.0], \
			[0.0,KInt0, -(KInt1 + self.perm1), Q1,0.0,0.0,0.0,0.0,0.0], \
			[0.0,0.0,KInt1, -(KInt2 + Q1 + self.perm2),0.0,0.0,0.0,0.0,0.0], \
			[0.0,0.0,0.0,KInt2, -(KInt3 + self.perm3), Q2,0.0,0.0,0.0], \
			[0.0,0.0,0.0,0.0,KInt3, -(KInt4 + Q2 + self.perm4),0.0,0.0,0.0], \
			[0.0,0.0,0.0,0.0,0.0,KInt4, -(KInt5 + self.perm5), Q3,0.0], \
			[0.0,0.0,0.0,0.0,0.0,0.0,KInt5, -(Kie + Q3 + self.perm6),0.0], \
			np.concatenate(([0.0],self.PermRates,[self.Kpel]))],'float')
		return Fmat

	def solveSys(self,t):
		self.t = t
		# Solve system of equations
		start_time = time.time()
		result = odeint(self.dF,self.M0,t,rtol=1e-3, atol=1e-3) #mxstep=500,hmax=100,Dfun=self.Fjac,full_output=True,
		#print("--- %s seconds ---" % np.str(time.time() - start_time))
		return result
	
	

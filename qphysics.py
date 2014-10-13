from visual import vector, mag, cross, mag, norm
import visual as vis
from functools import *
from math import *

DEBUG = False

# constants
#-------------------------
qe = 1.602176e-19 #C
me = 9.109382e-31 #kg
mp = 1.672621e-27 #kg
kc = 8.987552e9 #N m^2/C^2
G  = 6.673000e-11 #N m^2/kg^2 

rh = 5.3e-11 # m (H radius)

class Particle():

    def __init__(self,name,charge,mass):
        self.name = name
        self.charge = charge
        self.mass = mass

PARTICLES = {
    "e-":Particle("e-",-1.*qe,me),
    "e+":Particle("e+",+1.*qe,me),
    "p":Particle("p",qe,mp),
    "pbar":Particle("pbar",-1.*qe,mp)
    }

RSIZE = 0.4
LSIZE = 10.
MINR = 1e-6 # mm # minimum distance to be 'inside'

COLORQPOS = vis.color.blue # color positive charge
COLORQNO = vis.color.white # color no charge
COLORQNEG = vis.color.red # color negative charge

#--- utils

def sum_vectors(fs):
    f = reduce(lambda x,y:x+y,fs)
    return f

def efields(pos,qs):
    es = list(map(lambda q: q.efield(pos),qs))
    return es

def efield(pos,qs):
    es = efields(pos,qs)
    ee = sum_vectors(es)
    return ee

def set_color_charge(q):
    if (q.charge<0): q.vpy.color = COLORQNEG
    elif (q.charge>0): q.vpy.color = COLORQPOS
    else: q.vpy.color = COLORQNO

#--------- Qpys --------------------------------        

class Qpy:
    
    def __init__(self,vol,charge=me, mass=me, fix = False, mother = None,
                 velocity=vector(0,0,0), wrotation=vector(0,0,0)):
        self.vpy = vol
        self.time = 0.
        self.set_charge(charge)
        self.charge = charge
        self.mass = mass
        self.Imass = 0
        self.velocity = velocity
        self.wrotation = vector(0,0,0)
        self.fix = fix
        self.set_mother(mother)
        self.force = vector(0,0,0)
        self.torque = vector(0,0,0)
        return

    def set_charge(self,charge):
        self.charge = charge
        set_color_charge(self)

    def set_mother(self,mother):
        self.mother = mother
        if (mother == None): return
        self.fix = True
        if (DEBUG): print(id(self),' fix to ',id(mother))

    def is_inside(self,pos):
        ok = (mag(self.position()-pos)<MINR)
        if (DEBUG): print(id(self),' is inside ',pos,' ? ',ok)
        return ok

    def position(self):
        pos = self.vpy.pos
        if (self.mother): pos = self.mother.frame_to_world(pos)
        if (DEBUG): print(id(self),' position ',pos)
        return pos 

    def momentum(self):
        p = self.mass*self.velocity
        if (DEBUG): print(id(self),' momentum ',pos)
        return p

    def angular_momentum(self,ref=vector(0,0,0)):
        dir = self.position()-ref
        p = self.momentum()
        ll = cross(dir,p)
        if (DEBUG): print(id(self),' angular momentum ',ll)
        return ll

    def venergy(self):
        vv = mag(self.velocity)
        ee = 0.5*self.mass*vv*vv
        return ee

    def eenergy(self,qs):
        vv = sum(list(map(lambda q: q.vfield(self.position()),qs)))
        uu = self.charge*vv
        return uu        

    def external_forces(self,qps):
        # qs = list(filter(lambda q: id(q) != id(self) ,qps))
        es = list(map(lambda q: q.efield(self.position()),qps))
        fs = list(map(lambda ee: self.charge*ee,es))
        self.forces = fs
        if (DEBUG): print(id(self),' fs (N) ',fs)
        return fs

    def external_force(self,qps):
        fs = self.external_forces(qps)
        ff = sum_vectors(fs); self.force = ff
        if (DEBUG): print(id(self),' f (N) ',ff)
        return ff

    def external_torques(self,qps):
        # qs = list(filter(lambda q: id(q) != id(self) ,qps))
        ts = [vector(0,0,0)]; self.torques = ts
        if (DEBUG): print(id(self),' taus (Nm) ',ts)
        return ts

    def external_torque(self,qps):
        ts = self.external_torques(qps)
        t = sum_vectors(ts); self.torque = t
        if (DEBUG): print(id(self),' tau (Nm) ',t)
        return t
    
    def interaction(self,qps):
        qs = list(filter(lambda q: id(q) != id(self) ,qps))
        self.external_force(qs)
        self.external_torque(qs)
        if (DEBUG): 
            print(id(self),' F ',self.force,' (N),  tau ',self.torque,' (Nm) ')
        return 

    def move(self,dt):
        ac = (1./self.mass)*self.force
        self.vpy.pos = self.position() + dt*self.velocity
        self.velocity = self.velocity + dt*ac  
        if (DEBUG):
            print(id(self)," move F = ",self.force,"(N), a = ",ac
                  ," (m/s2), v = ",self.velocity,
                  " (m/s), x = ",self.position()," (m) ") 
        return 
            
    def turn(self,dt):
        vphi = self.wrotation*dt
        phi = mag(vphi); vaxis = vector(0,0,0)
        if (phi>0): 
            vaxis = vphi/phi 
            self.vpy.rotate(angle=phi,axis=vaxis)
        if (self.Imass==0): return
        ar = self.torque/self.Imass
        self.wrotation = self.wrotation + ar*dt
        if (DEBUG):
            print(id(self),' turn tau ',self.torque,' (Nm) ar ',ar,
                  ' (rad/s2) phi ',phi,' (rad) axis ',vaxis,
                  ' wrot ',self.wrotation,' (rad/s) ')
        return 
    
    def step(self,dt):
        self.time = self.time + dt
        if (DEBUG): print(id(self),' step ',dt,' (s) ')
        if (self.fix):
            if (DEBUG): print(id(self),' fix ',self.fix)
            return
        self.move(dt)
        self.turn(dt)
        return

class Charge(Qpy):

    def __init__(self,id=None, fix=False, 
                 charge=qe, mass=me, mother = None,
                 position=vector(0,0,0),  velocity=vector(0.,0.,0),
                 radius = RSIZE):
        sp = vis.sphere(pos=position, radius= radius, frame = mother)
        Qpy.__init__(self,sp, charge=charge, mass=mass, mother = mother, 
                     velocity = velocity)
        if (id != None):
            self.set_charge(PARTICLES[id].charge)
            self.mass = PARTICLES[id].mass
        self.fix = fix
        self.velocity = velocity

    def efield(self,r):
        if (self.is_inside(r)): return vector(0.,0.,0.)
        dir = r-self.position()
        dis = mag(dir)
        if (dis == 0.): return vector(0.,0.,0.)
        E = (kc*self.charge/(dis*dis))*(dir/dis)
        return E

    def vfield(self,r):
        if (self.is_inside(r)): return 0.
        dis = mag(r-self.position())
        if (dis == 0.): return 0.
        v = kc*self.charge/dis
        return v

class Wall(Qpy):
    
    def __init__(self,charge=1e6*qe,mass=1e6,size=vector(RSIZE,LSIZE,LSIZE),
                 fix=True, position=vector(0.,0.,0.),axis=vector(1.,0.,0.)):
        wl = vis.box(pos=position, axis=axis,
                     length=size.x,height=size.y,width=size.z)
        Qpy.__init__(self,wl,charge=charge,mass=mass,fix=fix)
        self.axis = norm(axis)
        self.size = size
        self.fix = fix
        area = self.size.y*self.size.z # size is along axis!
        self.sigma = self.charge/area # carefull!
        ee = mag(self.efield(vector(0,0,0)))
        print('Wall: charge ',self.charge, ' (C), area ',area, ' (m2)')
        print('Wall: sigma ',self.sigma, ' (C/m2), E ',ee, ' (N/C)')
        return

        #self.set_charge(charge)

    def efield(self,r):
        if (self.is_inside(r)): return vector(0.,0.,0.)
        E = 2.*pi*kc*self.sigma*self.axis
        # print('pi,kc,sigma,e',pi,kc,self.sigma,E,self.axis)
        return E        

    def vfield(self,r):
        if (self.is_inside(r)): return 0.
        V =  -2*pi*kc*self.sigma*(r-self.position()).dot(self.axis)
        return V

    def is_inside(self,r):
        pos = self.position()
        return mag(r-pos)<MINR
        x = abs(r.x-pos.x)
        y = abs(r.y-pos.y)
        z = abs(r.z-pos.z)
        #ok = ((x-self.size.x) <0 and (y-self.size.y)<0 and (z-self.size.z)<0)
        #return ok

class QComposite(Qpy):

    def __init__(self, vol, charge=me, mass=me, fix = False, 
                 velocity=vector(0,0,0), wrotation=vector(0,0,0)):
        Qpy.__init__(self, vol = vol, charge = charge, mass = me, fix = fix,
                     velocity = velocity, wrotation = wrotation)
        self.qs = []

    def add(self,qp):
        qp.set_mother(self.vpy)
        self.qs.append(qp)

    def is_inside(self,pos):
        oks = list(map(lambda q: q.is_inside(pos),self.qs))
        ok = sum(oks)>0
        return ok
    
    def efield(self,r):
        es = map(lambda q: q.efield(r),self.qs)
        ee = sum_vectors(es)
        if (DEBUG): print(id(self),' E ',ee, ' (N/C) ')
        return ee

    def vfield(self,r):
        vs = map(lambda q: q.vfield(r),self.qs)
        vv = sum(vs)
        if (DEBUG): print(id(self),' V ',vv, ' (V) ')
        return vv

    def external_forces(self,qs):
        fs = list(map(lambda q: q.external_forces(qs),self.qs))
        fs = list(reduce(lambda x,y:x+y,fs))
        if (DEBUG): print(id(self),' fs ',fs)
        return fs

    def external_torques(self,qs):
        fs = list(map(lambda q: q.external_force(qs),self.qs))
        def tau(q,f):
            dis = q.position()-self.position()
            t = cross(dis,f)
            if (DEBUG):
                print(id(self),'  t ',dis,f,t)
            return t
        ts = list(map(tau,self.qs,fs))
        if (DEBUG): print(id(self),' ts ',ts)
        return ts

class Dipole(QComposite):
    
    def __init__(self,distance=1.,charge=qe, mass = me, fix = False,
                 position=vector(0,0,0), axis=vector(1,0,0),
                 velocity = vector(0,0,0), wrotation = vector(0,0,0)):
        frame = vis.frame(pos=position,axis=axis)
        QComposite.__init__(self, frame, charge=0, mass=2.*mass, fix = fix,
                            velocity = velocity, wrotation = wrotation)
        pos = position+(distance/2.)*axis; xpos = frame.world_to_frame(pos)
        q1 = Charge(charge=charge,mass=mass,position=xpos,mother=frame)
        pos = position-(distance/2.)*axis; xpos = frame.world_to_frame(pos)
        q2 = Charge(charge=-1*charge,mass=mass,position=xpos,mother=frame)
        self.add(q1); self.add(q2) 
        xdir = vector(distance,0,0)
        self.tube = vis.cylinder(pos=xpos,axis=xdir,radius=RSIZE/2,
                                 color = vis.color.white, frame = frame)
        self.Imass = mass*(distance*distance)/2.

    def is_inside(self,r):
        ok1 = Qpy.is_inside(self,r)
        ok2 = QComposite.is_inside(self,r)
        ok = (ok1 or ok2)
        if (DEBUG): print(' is inside ',ok)
        return ok
                
    def external_forces(self,qps):
        qs =list(filter(lambda q: not q.is_inside(self.position()),qps))
        fs = QComposite.external_forces(self,qs)
        self.forces = fs
        return fs
    
    def external_torques(self,qps):
        qs =list(filter(lambda q: not q.is_inside(self.position()),qps))
        ts = QComposite.external_torques(self,qs)
        self.torques = ts
        return ts

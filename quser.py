from functools import *
from math import sqrt

import visual as vis
from visual.graph import gdisplay, gcurve
from visual import vector, mag, cross
import qpys.qphysics as qphys
import qpys.qvisual as qvis

DEBUG = False

# colors
COLORQPOS = vis.color.blue # color positive charge
COLORQNO = vis.color.white # color no charge
COLORQNEG = vis.color.red # color negative charge
COLORE = vis.color.magenta # color E-field
COLORF = vis.color.orange # color Force
COLORV = vis.color.yellow # color velocity

COLORG = vis.color.cyan # color velocity

MAXR = 0.002 # maximun distance to move in one step (m)
MAXPHI = 0.002 # maximum angle to turn in one step (rad)

#----------- Welcome
print(' Qpys: Electromagnetism with VPython')
print('\t\t J.A. Hernando, USC, 2013')

#----- Drawing elements

def axis(unit=1.,mode="xyz"):
    """ draw and return axis system.
    @ unit, indicate the size of the axis arrows
    @ mode = "xyz","yz","z",etc indicate what axis to draw 
    """
    xaxis  = []
    if (mode.find('x')>=0):
        xarrow = vis.arrow(pos=vector(0,0,0,),axis=vector(unit,0.,0.),
                       color=vis.color.green, shaftwidth = 0.2)
        xaxis.append(xarrow)
    if (mode.find('y')>=0):
        yarrow = vis.arrow(pos=vector(0,0,0,),axis=vector(0.,unit,0.),
                       color=vis.color.red, shaftwidth = 0.2)
        xaxis.append(yarrow)
    if (mode.find('z')>=0):
        zarrow = vis.arrow(pos=vector(0,0,0,),axis=vector(0.,0.,unit),
                       color=vis.color.blue, shaftwidth = 0.2)
        xaxis.append(zarrow)
    return xaxis

def grid(length=10.,unit=1.,mode="z",offset=0.):
    """ draw and return a grid in a plane
    @ length: the half width of the plane
    @ unit: the size of the squared to draw
    @ mode: "x","y","z" a plane in x, y, z
    @ offset: the offset of the plane of the grid (default in the origin 0.)
    """
    n = int(length/unit)
    xradius = 0.0
    lines = []
    for i in range(2*n+1):
        xline = vis.curve(color=vis.color.white,radius=xradius)
        yline = vis.curve(color=vis.color.white,radius=xradius)
        for j in range(2*n+1):
            x = (-length+i*unit,-length+j*unit,offset)
            y = (-length+j*unit,-length+i*unit,offset)
            if (mode == "x"): 
                x = (offset,-length+i*unit,-length+j*unit)
                y = (offset,-length+j*unit,-length+i*unit)
            elif (mode == "y"): 
                x = (-length+j*unit,offset,-length+i*unit)
                y = (-length+i*unit,offset,-length+j*unit)
            xline.append(pos=x)
            yline.append(pos=y)
        lines.append(xline); lines.append(yline)
    return lines

def surfpoints(f,unit=1,length=8):
    vs = []
    n = int(length/(2*unit))
    for i in range(n):
        for j in range(n):
            x0,y0 = -length/2+unit*i,-length/2+unit*j
            vs.append((x0,y0,f(x0,y0)))
    return vs

def maxsurfpoints(f,unit=1,length=8):
    vs = surfpoints(f,unit,length)
    v = max(list(map(lambda x: abs(x[2]),vs)))
    return v

def surfgrid(fun,unit=1,length=8,color=COLORG,scale=None):
    if (not scale):
        vs = maxsurfpoints(fun,unit,length)
        if (vv>0.): scale = length/(2.*vv)
    n = int(length/unit)
    xradius = 0.0
    lines = []
    for i in range(2*n+1):
        xline = vis.curve(color=color,radius=xradius)
        yline = vis.curve(color=color,radius=xradius)
        for j in range(2*n+1):
            x0,y0 = -length+i*unit,-length+j*unit
            x = (x0,y0,fun(x0,y0)*scale)
            x0,y0 = -length+j*unit,-length+i*unit
            y = (x0,y0,fun(x0,y0)*scale)
            xline.append(pos=x)
            yline.append(pos=y)            
        lines.append(xline); lines.append(yline)
    return lines

def draw_arrow(pos,dir,color):
    ar = vis.arrow(pos=pos,axis=dir,color=color,shaftwidth=0.2)
    return ar

def draw_cylinder(pos,dir,color,radius=0.1):
    ar = vis.cylinder(pos=pos,axis=dir,color=color,radius=radius)
    return ar

def draw_sphere(pos,color,radius=0.2):
    ar = vis.sphere(pos=pos,color=color,radius=radius)
    return ar

def color_charge(qob):
    if (qob.charge>0): return COLORQPOS
    elif (qob.charge<0): return COLORQNEG
    return COLORQNO

#---------------------------------------------
class Manager():
#---------------------------------------------

    def __init__(self,length,name="",unit=1):
        #self.scene = vis.display(title=name)
        #self.scene.select()
        if (name): vis.scene.title=name
        vis.scene.range = 1.2*length
        self.unit = unit
        self.length = length
        self.grid = grid(length,unit)
        self.maxr = MAXR # maximum step in distance allowed
        self.maxphi = MAXPHI # maximum angle allowed
        self.rate_slow = None # rate of program (large number (i.e 1000 = fast)
        self.rate_draw = 1
        self.qs = []
        self.draws = []
        self.store = []
        self.tup = []
        self.ftime = 1.
 
    def charge(self,id=None,charge=qphys.qe,mass=qphys.me, fix = False,
               position=vector(0,0,0), velocity = vector(0,0,0), 
               radius = qphys.RSIZE):
        q = qphys.Charge(id=id,charge=charge,mass=mass,fix=fix,
                         radius = radius,
                         position = position, velocity = velocity)
        self.add(q)
        return q

    def wall(self,charge=(1e6)*qphys.qe,mass=100, fix = True,
             position = vector(0,0,0), axis = vector(-1,0,0),
             size = vector(qphys.RSIZE,qphys.LSIZE,qphys.RSIZE)):
        w = qphys.Wall(charge=charge,mass=mass,fix=fix,
                       position=position,axis=axis,size=size)
        self.add(w)
        return w

    def dipole(self,distance=2.,charge=qphys.qe,mass=qphys.me, fix = False,
               position = vector(0,0,0), axis = vector (1,0,0,),
               velocity = vector(0,0,0), wrotation = vector(0,0,0)):
        d = qphys.Dipole(distance=distance,charge=charge,mass=mass,fix=fix,
                         position=position,axis=axis,
                         velocity=velocity,wrotation=wrotation)
        self.add(d)
        return d

    def add(self,q):
        self.qs.append(q); 
        #if (len(self.qs)>1): self.scales();
        self.scales();

    def qs_all(self):
        """ return all qs including the qs in composites 
        (i.e the charges of a dipole)
        """
        qs = self.qs
        for q in self.qs: 
            if (hasattr(q,'qs')): qs = qs + q.qs    
        return qs

    def qs_except(self,q0):
        """ return the qs except the argument (q0)
        """
        qs = []
        for q in self.qs: 
            if (id(q) != id(q0)): qs.append(q)
        return qs
    
    def qs_isin(self,q0):
        for q in self.qs: 
            if (id(q) == id(q0)): return True
        return False

    def scales(self):
        """ Compute the scale for the voltage, efield, force, velocity, time
        of a given configuration
        """
        L = 1.2*self.length
        self.uscale = 1; self.escale = 1;
        self.fscale = 1; self.vscale = 1; self.dt = 1;
        qs = self.qs_all()
        fun = lambda x,y: self.vfield(vector(x,y,0.))
        vv = maxsurfpoints(fun,self.unit,self.length)
        self.uscale = L/(2*vv)
        ee = max(list(map(lambda q: mag(qphys.efield(q.position(),qs)),qs)))
        if (ee>0): self.escale = L/(2.*ee)
        ff = max(list(map(lambda q: mag(q.external_force(qs)),qs)))
        if (ff>0): self.fscale = L/(2*ff)
        ac = max(list(map(lambda q: mag(q.external_force(qs))/q.mass,qs)))
        if (ac>0.): 
        #dtime = sqrt(2.*L/ac)/1000.
            self.dt = MAXR*sqrt(2.*L/ac)/10.
            self.vscale = L/(2.*sqrt(2.*L*ac))
        if (DEBUG):
            print('Vmax',vv,' Emax ',ee,' Fmax ',ff,' acmax ',ac,
                  ' vmax ',L/(2*vscale))
            print('vscale',self.uscale,' (V) escale ',self.escale,
                  ' (N), fscale ',self.fscale,
                  ' (N/m) vscale ',self.vscale,' (1/s) dt ',
                  self.dt,' (s) ')
        return

    def input(self):
        qs = dialogue(self)
        return qs

    def efield(self,pos):
        ee = qphys.efield(pos,self.qs)
        return ee

    def efields(self,pos):
        es = qphys.efields(pos,self.qs)
        return es

    def vfield(self,pos):
        vv = sum(list(map(lambda q: q.vfield(pos),self.qs)))
        return vv        

    def eenergy(self,q0):
        return q0.energy(selg.qs)

    def draw(self):
        for dr in self.draws: dr()

    def draw_forces(self):
        ars = []
        for qi in self.qs:
            ars.append(self.draw_forces(qi))
        return ars

    def draw_forces(self,qob,step=True):
        qob.interaction(self.qs)
        pos = qob.position(); fs = qob.forces
        ars = list(map(lambda f : draw_arrow(pos,self.fscale*f,COLORF),fs))
        qob.arrow_forces = ars
        def draw():
            for i in range(len(qob.forces)):
                qob.arrow_forces[i].pos = qob.position()
                qob.arrow_forces[i].axis = qob.forces[i] * self.fscale
            return
        draw()
        if (step): self.draws.append(draw) 
        return ars

    def draw_force(self,qob,step=True):
        qob.interaction(self.qs);
        ar = draw_arrow(qob.position(),qob.force*self.fscale,COLORF)
        qob.arrow_force = ar
        def draw():
            qob.arrow_force.pos = qob.position()
            qob.arrow_force.axis = qob.force * self.fscale
        draw();
        if (step): self.draws.append(draw)
        return ar

    def draw_velocity(self,qob,step=True):
        ar = draw_arrow(qob.position(),qob.velocity*self.vscale,COLORV)
        qob.arrow_velocity = ar
        def draw():
            qob.arrow_velocity.pos = qob.position()
            qob.arrow_velocity.axis = qob.velocity*self.vscale
        draw()
        if (step): self.draws.append(draw)
        return ar

    def draw_trajectory(self,qob):
        trj = vis.curve(color=qob.vpy.color)
        qob.trajectory = trj
        def draw():
            qob.trajectory.append(pos=qob.position())
        draw()
        self.draws.append(draw)
        return trj

    def draw_torques(self,qob):
        qob.interaction(self.qs)
        pos = qob.position(); fs = qob.torques
        ars = list(map(lambda f : draw_arrow(pos,self.fscale*f/(0.5*self.length),COLORF),fs))
        qob.arrow_torques = ars
        def draw():
            for i in range(len(qob.forces)):
                qob.arrow_torques[i].pos = qob.position()
                qob.arrow_torques[i].axis = qob.torques[i] * self.fscale/(0.5*self.length)
            return
        draw()
        if (step): self.draws.append(draw) 
        return ars
        
    def draw_torque(self,qob,step=True):
        qob.interaction(self.qs);
        ar = draw_arrow(qob.position(),qob.torque*self.fscale/(0.5*self.length),COLORF)
        qob.arrow_torque = ar
        def draw():
            qob.arrow_torque.pos = qob.position()
            qob.arrow_torque.axis = qob.torque * self.fscale/(0.5*self.length)
        draw();
        if (step): self.draws.append(draw)
        return ar

    def draw_wrotation(self,qob,step=True):
        ar = draw_arrow(qob.position(),qob.wrotation*self.vscale*self.length/2.,COLORV)
        qob.arrow_wrotation = ar
        def draw():
            qob.arrow_wrotation.pos = qob.position()
            qob.arrow_wrotation.axis = qob.wrotation*self.vscale*self.length/2.
        draw()
        if (step): self.draws.append(draw)
        return ar

    def draw_energy(self,qob,step=True):
        scale = self.uscale/abs(qob.charge)
        qs = self.qs_except(qob)
        uv = qob.venergy()*scale; ue = qob.eenergy(qs)*scale
        color = color_charge(qob)
        are = draw_cylinder(qob.position(),ue*vector(0,0,1),color)
        qob.cylinder_ee = are
        arv = draw_cylinder(qob.position()+ue*vector(0,0,1),
                            uv*vector(0,0,1),COLORV,radius=0.12)
        qob.cylinder_ev = arv
        arq = draw_sphere(qob.position()+ue*vector(0,0,1),color,radius=0.2)
        qob.sphere_ene = arq
        def draw():
            uv = qob.venergy()*scale; ue = qob.eenergy(qs)*scale
            qob.cylinder_ee.pos = qob.position()
            qob.cylinder_ee.axis = ue*vector(0,0,1)
            qob.cylinder_ev.pos = qob.position()+ue*vector(0,0,1)
            qob.cylinder_ev.axis = uv*vector(0,0,1)
            qob.sphere_ene.pos = qob.position()+ue*vector(0,0,1)
        draw()
        if (step): self.draws.append(draw)
        return are,arv,arq

    def draw_efields(self,pos):
        es = qphys.efields(pos,self.qs)
        ars = list(map(lambda ee: draw_arrow(pos,self.escale*ee,COLORE),es))
        for ar in ars: self.store.append(ar)
        return ars
        
    def draw_efield(self,pos):
        ee = qphys.efield(pos,self.qs)
        ar = draw_arrow(pos,self.escale*ee,COLORE)
        self.store.append(ar)
        return ar

    def draw_efield_at(self,qob,step=True):
        pos = qob.position(); ee = self.efield(pos)
        ar = draw_arrow(pos,ee*self.escale,COLORE)
        qob.arrow_efield = ar
        def draw():
            pos = qob.position(); ee = self.efield(pos)
            qob.arrow_efield.pos = pos
            qob.arrow_efield.axis = ee*self.escale
        draw()
        if (step): self.draws.append(draw)
        return ar

    def draw_vfield(self,unit=1,factor=1.):
        f = lambda x,y: factor*self.vfield(vector(x,y,0.))
        color = COLORQPOS
        if (factor<0): color = COLORQNEG
        ss = surfgrid(f,unit,self.length,color=color,scale=self.uscale)
        self.store.append(ss)
        return ss


    def draw_vfield_on(self,q0,unit=1):
        qs = self.qs_except(q0)
        factor = 1.
        if q0.charge<0: factor = -1.
        f = lambda x,y: factor*sum(list(map(lambda q: q.vfield(vector(x,y,0)),qs)))
        color = color_charge(q0)
        ss = surfgrid(f,unit,self.length,color=color,scale=self.uscale)
        self.store.append(ss)
        return ss

    def clear(self):
        for vi in self.store : vi.visible = False

    def delete(self):
        for q in self.qs: q.visible = False
        self.qs = []

    def interaction(self):
        for q in self.qs: q.interaction(self.qs)
    
    def dtime(self,dt=-1.):
        dt0 = max(dt,self.dt)
        fs1 = list(map(lambda q: mag(q.velocity*dt0)/self.maxr,self.qs))
        fs2 = list(map(lambda q: mag(q.wrotation*dt0)/self.maxphi,self.qs))
        fs = fs1+fs2; f = max(fs)
        if (f>1.): 
            self.ftime = self.ftime/f
            dt = dt0/f
            self.dt = dt
        else: dt = dt0
        if (DEBUG): 
            print(' distances ',list(map(lambda x: x*self.maxr,fs1)))
            print(' phis ',list(map(lambda x: x*self.maxphi,fs2)))
            print(' dtime ',dt0,' -> ',dt)
        return dt

    def step(self,dt=-1):
        """ makes an step in time (dt), all the qobs are involved in the step
        """
        dt = self.dtime(dt)
        list(map(lambda qo: qo.step(dt),self.qs))
        return

    def run(self,nsteps=1,dtime=-1,fun=None,rate_fun=1):
        """ run for a time the interaction between all qobs,
        each step has a dtime (delta time), after every step the function
        fun is called with qobs as arguments, 
        every rate events the display is refresh
        """
        for istep in range(nsteps): 
            if (self.rate_slow): vis.rate(self.rate_slow)
            self.interaction(); self.step(dtime); 
            if (istep%self.rate_draw == 0): 
                for dr in self.draws: dr()
            if ((fun != None) and (istep%rate_fun == 0)): 
                self.tup.append(fun())
        if (self.ftime<1): print('Set time factor ',self.ftime)
        return self.tup
    
    def xydrag_follow(self,q,fun=None):
        self.tup = []
        def xyfollow():
            q.vpy.pos.z = 0 # problem with z y drawing... 
            self.interaction() # interaction
            self.draw() # draw 
            if (fun):  
                val = fun() # execute function
                #print('ret ',val)
                self.tup.append(val) # store in tuple    
        qvis.drag_object(q.vpy,xyfollow)
        return self.tup

    def xyclick(self,fun=None):
        vals = []
        def xyfun(xpos):
            xpos.z = 0.
            if (abs(xpos.x)>self.length or abs(xpos.y)>self.length): 
                return False
            if (fun):
                val = fun(xpos); vals.append(val)
            return True
        ok = True
        while ok:
            ok = qvis.click_position(xyfun)
        return vals

    def xyclick_draw_efield(self):
        return self.xyclick(self.draw_efield)
    
    def xyclick_draw_efields(self):
        return self.xyclick(self.draw_efields)
    
#----------------------------------
#    Interaction with the Client 
#----------------------------------

def dialogue(man):
    qs = []
    ok = True
    while (ok):
        comment = ' Input 0-exit, 1-charge, 2-wall, 3-dipole : '
        val = int(input(comment))
        if (val <=0): ok = False
        elif (val == 1):
            comment = ' Enter: x,y,z,q,t position, charge -in e units- and t = 0/1 for free/fix charge: '
            val = input(comment)
            print('your input ',val)
            vals = val.split(',')
            x,y,z,c,t = list(map(float,vals))
            print('charge ',c)
            q = man.charge(position=vector(x,y,z),
                           charge=c*qphys.qe,mass=qphys.me,fix=bool(t))
            qs.append(q)        
    return qs

#--------- delaing with graphs

def fun_retrieve(q,vals):
    funs = []
    def fun(val):
        if (val == "t"): return lambda q: q.time
        if (val == "x"): return lambda q: q.position().x
        if (val == "y"): return lambda q: q.position().y
        if (val == "z"): return lambda q: q.position().z
        if (val == "dd"): return lambda q: mag(q.position())
        if (val == "fx"): return lambda q: q.force.x
        if (val == "fy"): return lambda q: q.force.y
        if (val == "fz"): return lambda q: q.force.z
        if (val == "ff"): return lambda q: mag(q.force)
        if (val == "vx"): return lambda q: q.velocity.x
        if (val == "vy"): return lambda q: q.velocity.y
        if (val == "vz"): return lambda q: q.velocity.z
        if (val == "vv"): return lambda q: mag(q.velocity)
        if (val == "tx"): return lambda q: q.torque.x
        if (val == "ty"): return lambda q: q.torque.y
        if (val == "tz"): return lambda q: q.torque.z
        if (val == "tt"): return lambda q: mag(q.torque)
        if (val == "wx"): return lambda q: q.wrotation.x
        if (val == "wy"): return lambda q: q.wrotation.y
        if (val == "wz"): return lambda q: q.wrotation.z
        if (val == "ww"): return lambda q: mag(q.wrotation)
    for val in vals:
        funs.append(fun(val))
    def ret():
        zvals = list(map(lambda f: f(q),funs))
        return zvals
    return ret
        
def get_column(ntup,i):
    x = list(map(lambda t: t[i],ntup))
    return x

def plot(x,y,xtitle="",ytitle=""):
    gd = gdisplay(xtitle= xtitle,ytitle = ytitle)
    gf = gcurve(color = COLORG)
    for i in range(len(x)): gf.plot(pos=(x[i],y[i]))
    return [gd,gf]


def multiplot(xs,ys,xtitle="",ytitle=""):
    gd = gdisplay(xtitle= xtitle,ytitle = ytitle)
    gfs = []
    for k in range(len(ys)):
        k0 = -1
        if (k<len(xs)): k0 = k
        xx = xs[k0]; yy = ys[k]
        gf = gcurve(color = COLORG)
        for i in range(len(xx)): gf.plot(pos=(xx[i],yy[i]))
        gfs.append(gf)
    return [gd,gfs]

#------- mouse


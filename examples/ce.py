#
# Examples using VPython about the Electric Field (Chapter -I)
#
# Jose A. Hernando 18/12/2012, UCC
#

from visual.graph import gdisplay, gcurve
from visual import vector, mag, color
import qpys
from qpys import vector, mag
import visual as vis
from functools import *
from math import *

COLORG = color.cyan

def coulomb0():
    """
    """

    # create the grid
    L = 8 # size of the grid
    man = qpys.Manager(L,'Coulomb')

    # create the charges
    q1 = man.charge(id="e+",position=vector(0.,0.,0),fix=True)  
    q2 = man.charge(id="e+",position=vector(+1.,0.,0),fix=True)  
    q3 = man.charge(id="e+",position=vector(+2.,0.,0),fix=True)  
    q4 = man.charge(id="e+",position=vector(+3.,0.,0),fix=True)  
    
    # draw the forces on the charges
    man.draw_force(q1); man.draw_force(q2)
    print('force q2 ',q2.force)
    
    input('enter key')
    q2.vpy.position = vector(+2.,0.,0.)
    man.draw_force(q1); man.draw_force(q2)
    print('force q2 ',q2.force)

    return man


#--------------------------------------------------
def coulomb():
#--------------------------------------------------
    """ Example: Coulomb charge
    Creates two charges, one positive fix (positron) and a second one that moves
    It ask for trhee time Ask for a charge for the second charge (in units of electron). 
    After you enter the charge you can  drag the charge to show the Coulomb forces, click anywere in the panel to continue.
    After tree trials, a graph with the force vs distance for the 3 charges that you enter.
    @ author: J.A. Hernando , date 18/12/2012
    """

    # create the grid
    L = 8 # size of the grid
    man = qpys.Manager(L,'Coulomb')

    # create the charges
    q1 = man.charge(id="e+",position=vector(0.,0.,0),fix=True)  
    q2 = man.charge(id="e+",position=vector(+4.,0.,0),fix=True)  

    # draw the forces on the charges
    man.draw_force(q1); man.draw_force(q2)

    # get a function to retrieve distance and force at each step
    ret = qpys.fun_retrieve(q2,['dd','ff']) 
    print('Drag and move the charge around. Click elsewhere to exit.')
    # drag the 2nd charge, move around, distance and force are stores into a ntuple
    tup = man.xydrag_follow(q2,ret) 

    # get the values and plot them
    dd = qpys.get_column(tup,0)
    ff = qpys.get_column(tup,1)
    cc = qpys.plot(dd,ff,"distance (m) ","force (N)")

    print("That's all folks!")
    return man

#--------------------------------------------------
def coulomb2():
#--------------------------------------------------
    """ Example: Coulomb charge
    Creates two charges, one positive fix (positron) and a second one that moves
    It ask for trhee time Ask for a charge for the second charge (in units of electron). 
    After you enter the charge you can  drag the charge to show the Coulomb forces, click anywere in the panel to continue.
    After tree trials, a graph with the force vs distance for the 3 charges that you enter.
    @ author: J.A. Hernando , date 18/12/2012
    """

    # create the grid
    L = 8 # size of the grid
    man = qpys.Manager(L,'Coulomb')

    # create the charges
    q1 = man.charge(id="e+",position=vector(0.,0.,0),fix=True)  
    ntrials = 2
    dds = []; ffs = []
    q2 = man.charge(id='e+',position=vector(+4.,0.,0),fix=True)  
    # draw the forces on the charges
    man.draw_force(q1); man.draw_force(q2)
    for i in range(ntrials):
        if (i>0):
            qq = float(input(' Enter charge (in e units) '))
            q2.vpy.radius = qpys.RSIZE*sqrt(abs(qq))
            q2.set_charge(qq*qpys.qe)
        man.interaction()
        man.draw()

        # get a function to retrieve distance and force at each step
        ret = qpys.fun_retrieve(q2,['dd','ff']) 
        print('Drag and move the charge around. Click elsewhere to exit.')
        # drag the 2nd charge, move around, distance and force are stores into a ntuple
        tup = man.xydrag_follow(q2,ret) 
        dds.append(qpys.get_column(tup,0))
        ffs.append(qpys.get_column(tup,1))

    # plot the values
    gfs = qpys.multiplot(dds,ffs,' distance (m)',' force (N)')

    print("That's all folks!")
    return man

#--------------------------------------------------------
def sumforces1():
#--------------------------------------------------------
    """ Example of adding forces of different charges.
    Three positive charges in a line. What happend to the one in the middle?
    """
    L = 8
    man = qpys.Manager(L,'Trap')

    q1 = man.charge(id="e+",position=vector(-6,0,0),fix=True)
    q2 = man.charge(id="e+",position=vector(+6,0,0),fix=True)
    q3 = man.charge(id="e+",position=vector(3,0,0),fix = False)
    qs = [q1,q2,q3]

    man.fscale = 1.4 * man.fscale
    print('All charges are equal')
    input("What are the forces acting on the charges? Enter key ")
    ars = list(map(lambda q: man.draw_forces(q),qs))

    input("What are the equivalent force? Enter key ")
    for iars in ars:
        for ar in iars: ar.visible = 0
    ars = list(map(lambda q: man.draw_force(q),qs))
    
    #vscale = 0.5 # scale the velocity by this value
    comment = "The second charge can move. How? Enter key."
    input(comment)
    man.draw_velocity(q3)

    ret = qpys.fun_retrieve(q3,['t','fx','vx','x'])                
    man.rate_slow = 1000 # sloww the movement...
    tup = man.run(nsteps=20000,fun=ret,rate_fun=50)

    #print(len(tup),tup[0])
    
    ts = qpys.get_column(tup,0)
    fs = qpys.get_column(tup,1)
    vs = qpys.get_column(tup,2)
    xs = qpys.get_column(tup,3)

    c1 = qpys.plot(ts,fs,'time (s)','force (N)')
    c2 = qpys.plot(ts,vs,'time (s)','x-velocity (m/s)')
    c3 = qpys.plot(ts,xs,'time (s)','x-position (m)')
    print("That's all folks!")
    return man

#--------------------------------------------------------
def sumforces2():
#--------------------------------------------------------
    """ Example of adding forces of different charges.
    Three charges, two positives and one negative in the axis of the two.
    """

    L = 8
    man = qpys.Manager(L,'Trap')

    q1 = man.charge(id="e+",position=vector(-6,0,0),fix=True)
    q2 = man.charge(id="e+",position=vector(+6,0,0),fix=True)
    q3 = man.charge(id="e-",position=vector(0,6,0),fix = False)
    qs = [q1,q2,q3]

    print('All charges are equal.')
    input("What are the forces acting on the charges? Enter key.")
    ars = list(map(lambda q: man.draw_forces(q),qs))

    input("What are the equivalent force? Enter key ")
    for iars in ars:
        for ar in iars: ar.visible = 0
    ars = list(map(lambda q: man.draw_force(q),qs))
    
    #vscale = 0.5 # scale the velocity by this value
    comment = "Negative charge will move. How? Enter key."
    input(comment)
    man.draw_velocity(q3)

    ret = qpys.fun_retrieve(q3,['t','fy','vy','y'])                
    man.rate_slow = 1000 # sloww the movement...
    tup = man.run(nsteps=20000,fun=ret,rate_fun=50)

    ts = qpys.get_column(tup,0)
    fs = qpys.get_column(tup,1)
    vs = qpys.get_column(tup,2)
    xs = qpys.get_column(tup,3)

    c1 = qpys.plot(ts,fs,'time (s)','force (N)')
    c2 = qpys.plot(ts,vs,'time (s)','y-velocity (m/s)')
    c3 = qpys.plot(ts,xs,'time (s)','y-position (m)')
    print("That's all folks!")
    return man

#--------------------------------------------------------
def sumforces3():
#--------------------------------------------------------
    """ Example of adding forces of different charges.
    For equal charges in a square. One will be free and move.
    """
    L = 8
    man = qpys.Manager(L,'Trap')

    d = 5
    q1 = man.charge(id="e+",position=vector(0,d,0),fix=True)
    q2 = man.charge(id="e+",position=vector(0,-1.*d,0),fix=True)
    q3 = man.charge(id="e-",position=vector(d,0,0),fix = True)
    q4 = man.charge(id="e-",position=vector(-1.*d,0,0),fix = False)
    
    qs = [q1,q2,q3,q4]
    man.fscale = 0.5*man.fscale

    print('All charges as equal')
    input("What are the forces acting on the charges? Enter key.")
    ars = list(map(lambda q: man.draw_forces(q),qs))

    input("What are the equivalent force? Enter key ")
    for iars in ars:
        for ar in iars: ar.visible = 0
    ars = list(map(lambda q: man.draw_force(q),qs))
    
    comment = "The negative charge on the left will move. How? Enter key"
    input(comment)
    man.draw_velocity(q4)

    ret = qpys.fun_retrieve(q4,['t','fx','vx','x'])                
    man.rate_slow = 1000 # sloww the movement...
    tup = man.run(nsteps=20000,fun=ret,rate_fun=50)

    ts = qpys.get_column(tup,0)
    fs = qpys.get_column(tup,1)
    vs = qpys.get_column(tup,2)
    xs = qpys.get_column(tup,3)

    c1 = qpys.plot(ts,fs,'time (s)','force (N)')
    c2 = qpys.plot(ts,vs,'time (s)','velocity (m/s)')
    c3 = qpys.plot(ts,xs,'time (s)','position (m)')
    print("That's all folks!")
    return man

#--------------------------------------------------------
def sumforces4():
#--------------------------------------------------------
    """ Example of adding forces of different charges.
    For equal charges in a square. One will be free and move.
    """
    L = 8
    man = qpys.Manager(L,'Trap')

    d = 5
    q1 = man.charge(id="e+",position=vector(0,d,0),fix=True)
    q2 = man.charge(id="e+",position=vector(0,-1.*d,0),fix=True)
    q3 = man.charge(id="e-",position=vector(d,0,0),fix = False)
    q4 = man.charge(id="e-",position=vector(-1.*d,0,0),fix = False)
    
    qs = [q1,q2,q3,q4]

    ars = list(map(lambda q: man.draw_force(q),qs))
    
    comment = "The negative charges will move. How? Enter key"
    input(comment)
    man.draw_velocity(q3); man.draw_velocity(q4)

    man.rate_slow = 1000 # sloww the movement...
    tup = man.run(nsteps=30000)

    print("That's all folks!")
    return man

#----------------------------------------------------
def efield1():
#----------------------------------------------------
    """ Efield fron a electric charge
    """
    L = 8
    man = qpys.Manager(L,'Electric Field')
    q1 = man.charge('e+',position=vector(0,0,0),fix=True)
    #qs = man.input()
    
    # set Efield-scale by hand 
    man.escale = 1./(4e-11)
    
    print('Click in a position on the grid to see the electric field')
    print('Click outside the grid to exit.')
    vars = man.xyclick_draw_efield() 

    q2 = man.charge('e-',position=vector(0,5,0),fix=True)
    print('Drag the negative charge and move it around to see the electric field and the force acting on it')
    man.clear()

    man.draw_force(q2)
    man.draw_efield_at(q2)

    man.xydrag_follow(q2)
    print("That's all folks!") 
    return man 

#----------------------------------------------------
def efield2():
#----------------------------------------------------
    """ Efield from 2 positive chages
    """
    L = 8
    man = qpys.Manager(L,'Dipole Electric Field')
    q1 = man.charge('e+',position=vector(-3,0,0),fix=True)
    q2 = man.charge('e+',position=vector(3,0,0),fix=True)
    #qs = man.input()
        
    print('Click in a position on the grid to see the electric fields')
    print('Click outside the grid to exit.')
    vars = man.xyclick_draw_efields() 

    print('Click in a position on the grid to see the total electric field')
    print('Click outside the grid to exit.')
    man.clear()
    vars = man.xyclick_draw_efield() 

    q3 = man.charge('e-',position=vector(0,5,0),fix=True)
    print('Drag the negative charge and move it around to see the electric field and the force acting on it')
    man.clear()
    man.draw_force(q3)
    man.draw_efield_at(q3)
    man.xydrag_follow(q3)
    
    
    #print('Graph of the x-component of the electric field along the x coordinate')
    #n= 50; xstep = 2.*L/(1.*n)
    #xs = list(map(lambda i: xstep*i-L,range(n)))
    #es = list(map(lambda x: man.efield(vector(x,0,0)).x,xs))
    #c = qpys.plot(xs,es,' x (m)',' Ex (N/C)')

    print("That's all folks!") 
    return man 


#----------------------------------------------------
def efield3():
#----------------------------------------------------
    """ Efield from a dipole
    """
    L = 8
    man = qpys.Manager(L,'Dipole Electric Field')
    q1 = man.charge('e+',position=vector(-3,0,0),fix=True)
    q2 = man.charge('e-',position=vector(3,0,0),fix=True)
    #qs = man.input()
        
    print('Click in a position on the grid to see the electric fields')
    print('Click outside the grid to exit.')
    vars = man.xyclick_draw_efields() 
    man.clear()

    man.clear()
    print('Click in a position on the grid to see the total electric field')
    print('Click outside the grid to exit.')
    vars = man.xyclick_draw_efield() 
    man.clear()

    q3 = man.charge('e-',position=vector(0,5,0),fix=True)
    print('Drag the negative charge and move it around to see the electric field and the force acting on it')

    man.draw_force(q3)
    man.draw_efield_at(q3)
    man.xydrag_follow(q3)

    print("That's all folks!") 
    return man 
    
#-----------------------------------------------------------
def move1():
#-----------------------------------------------------------
    """ Sets few charges in a constant electic field.
    Show the Electric field and the movement of the charges.
    """
    L = 8
    man = qpys.Manager(L,'Infinite Wall Electric Field ')
    
    q1 = man.charge('e-',position=vector(-6,2,0,),fix=False)
    q2 = man.charge('e-',position=vector(-4,0,0,),fix=False)
    q3 = man.charge('e-',position=vector(-2,-2,0,),fix=False)
    
    w1 = man.wall(charge=1e-6, mass = 100, fix=True, 
                  position=vector(8,0,0),axis=vector(-1,0,0),
                  size=vector(0.4,16.,16.))

    print('Click to see the E field.')
    print('Click outiside the grid to continue.')
    man.xyclick_draw_efield()
    man.clear()

    def bounce():
        for q in [q1,q2,q3]:
            if (q.position().x >= w1.position().x):
                q.velocity = -1.*q.velocity
            
    input('How the charges will move? Enter key')
    man.dt = man.dt*0.45
    vals = man.run(nsteps=35000,fun=bounce,rate_fun=1)
    print("That's all folks!")
    return man

#----------------------------------------------------------------------
def move2():
#----------------------------------------------------------------------
    """ Sets few charges in a constant electic field.
    Show the Electric field and the movement of the charges.
    """
    L = 8
    man = qpys.Manager(L,'Infinity Wall Electric Field')
    
    q1 = man.charge('e-',fix=False,
                    position=vector(-6,0.5,0),velocity=vector(2e4,0,0))
    q2 = man.charge('e+',fix=False,
                    position=vector(-6,-0.5,0),velocity=vector(2e4,0,0))

    w1 = man.wall(charge=(1e6)*qpys.qe,mass = 100,
                  fix=True,size=vector(0.4,16.,16.),
                  position=vector(0,-8,0),axis=vector(0,1,0))

    qs = [q1,q2,w1]
    man.draw_trajectory(q1); man.draw_trajectory(q2)
    man.draw_velocity(q1); man.draw_velocity(q2)

    print('Click to see the E field.')
    print('Click outiside the grid to continue.')
    man.xyclick_draw_efield()
    man.clear()

    ee = man.efield(vector(0,0,0))
    man.draw_efield(vector(0,0,0))
    print('Charges are an electron and a positron.')
    print('Electric field (N/C) ',ee)
    print('Position e-: ',q1.position(),', e+: ',q2.position(), ' (m) ')
    print('Velocity e-: ',q1.velocity,', e+: ',q2.velocity, ' (m/s)')


    input('How the charges will move? Enter key.')
    man.dt = man.dt*0.6
    man.run(nsteps=8000)
    print('Position e-: ',q1.position(),', e+: ',q2.position(), ' (m) ')
    print('Velocity e-: ',q1.velocity,', e+: ',q2.velocity, ' (m/s) ')
    print("That's all folks!")
    return man
    
#-----------------------------------------------------------
def move3():
#-----------------------------------------------------------
    """ Sets few charges in a constant electic field.
    Set the value of the charge for both to reach at the same time to the plane
    """
    L = 8
    man = qpys.Manager(L,'Infinite Wall Electric Field')
    
    q1 = man.charge('e-',position=vector(-4,2,0,),fix=False)
    q2 = man.charge('e-',position=vector(2,-2,0,),fix=False)
    
    w1 = man.wall(charge=(1e6)*qpys.qe, mass = 100, fix=True, 
                  position=vector(8,0,0),axis=vector(-1,0,0),
                  size=vector(0.4,16.,16.))

    man.draw_efield(vector(0,0,0))
    #ee = man.efield(vector(0,0,0))
    #print('Efield ',ee, ' (N/C)')

    q = input('What must be the charge value of charge far away from the wall to arrive at the same time? Enter value (in e units): ')
    q1.set_charge(int(q)*qpys.qe)
                             
    def bounce():
        for q in [q1,q2]:
            if (q.position().x >= w1.position().x):
                q.velocity = -1.*q.velocity            
                
    ret1 = qpys.fun_retrieve(q1,['t','x'])
    ret2 = qpys.fun_retrieve(q2,['t','x'])
    def ret(): 
        bounce()
        return ret1()+ret2() 

    man.dt = man.dt*0.35
    vals = man.run(nsteps=60000,fun=ret,rate_fun=1)
    
    t1 = qpys.get_column(vals,0)
    x1 = qpys.get_column(vals,1)
    t2 = qpys.get_column(vals,2)
    x2 = qpys.get_column(vals,3)
    
    c1 = qpys.plot(t1,x1)
    c1 = qpys.plot(t2,x2)

    print("That's all folks!")
    return man

#-------------------------------------------------
def orbit1():
    """ Sets few charges in a constant electic field.
    Set the value of the charge for both to reach at the same time to the plane
    """
    L = 8
    man = qpys.Manager(L,'Kepler')
    
    f = qpys.me/qpys.mp # ratio of masses e-/proton
    q1 = man.charge('e-',position=vector(0,6,0,),velocity=vector(6.5,0,0),fix=False)
    q2 = man.charge('p',position=vector(0,0,0,),velocity=f*vector(-6.5,0,0),fix=False)

    man.draw_velocity(q1); man.draw_velocity(q2);
    print('System of an electron and a proton')
    print('e- velocity ',q1.velocity,' (m/s) ')
    input('How are the forces? Enter key')
    man.fscale = man.fscale/3.
    man.draw_force(q1); man.draw_force(q2)

    ret = qpys.fun_retrieve(q1,['t','dd','vv','ff'])

    input('How they will move? Enter key')
    man.dt = man.dt*0.90
    man.rate_slow = 1000
    man.draw_trajectory(q1); man.draw_trajectory(q2)
    tup = man.run(30000,fun=ret,rate_fun=300)

    print('len tup ',len(tup),tup[0],tup[-1])
    ts = qpys.get_column(tup,0)
    dd = qpys.get_column(tup,1)
    vv = qpys.get_column(tup,2)
    ff = qpys.get_column(tup,3)

    c1 = qpys.plot(ts,dd,'time (s)','distance (m)')
    c2 = qpys.plot(ts,vv,'time (s)',' velocity (m/s)')
    c3 = qpys.plot(ts,ff,'time (s)','force (N)')

    print("That's all folks!")
    return man

#-------------------------------------------------
def orbit2():
    """ Sets few charges in a constant electic field.
    Set the value of the charge for both to reach at the same time to the plane
    """
    L = 8
    man = qpys.Manager(L,'Kepler')
    
    f = qpys.me/qpys.mp # ratio of masses e-/proton
    q1 = man.charge('e-',position=vector(0,7,0,),velocity=vector(4.,0,0),fix=False)
    q2 = man.charge('p',position=vector(0,0,0,),velocity=f*vector(-4.,0,0),fix=False)

    man.draw_velocity(q1); man.draw_velocity(q2);
    print('System of an electron and a proton')
    print('e- velocity ',q1.velocity,' (m/s) ')
    input('How are the forces? Enter key')
    man.fscale = man.fscale/5.
    man.draw_force(q1); man.draw_force(q2)

    ret = qpys.fun_retrieve(q1,['t','x','y','vx','vy'])

    input('How they will move? Enter key')
    man.dt = man.dt*0.35
    man.rate_slow = 1000
    man.draw_trajectory(q1); man.draw_trajectory(q2)
    
   #tup = man.run(25000,fun=ret,rate_fun=300)
    tup = man.run(55000,fun=ret,rate_fun=300)

    ts = qpys.get_column(tup,0)
    xx = qpys.get_column(tup,1)
    yy = qpys.get_column(tup,2)
    vx = qpys.get_column(tup,3)
    vy = qpys.get_column(tup,4)

    c1 = qpys.plot(ts,xx,'time (s)','x distance (m)')
    c2 = qpys.plot(ts,yy,'time (s)','y distance (m)')
    c3 = qpys.plot(ts,vx,'time (s)',' vx velocity (m/s)')
    c4 = qpys.plot(ts,vy,'time (s)','vy velocity (m/s)')

    print("That's all folks!")
    return man,[c1,c2,c3,c4]

#-------------------------------------------------
def orbit3():
    """ Sets few charges in a constant electic field.
    Set the value of the charge for both to reach at the same time to the plane
    """
    L = 8
    man = qpys.Manager(L,'Kepler')
    
    q1 = man.charge('e-',position=vector(0,6,0,),velocity=vector(2,0,0),fix=False)
    q2 = man.charge('e+',position=vector(0,-6,0,),velocity=vector(-2,0,0),fix=False)

    man.draw_velocity(q1); man.draw_velocity(q2);
    print('System of an electron and a positron')
    print('e- velocity ',q1.velocity,' (m/s) ')
    input('How are the forces? Enter key')
    man.fscale = man.fscale/4.
    man.draw_force(q1); man.draw_force(q2)

    input('How they will move? Enter key')
    man.dt = man.dt*0.35
    man.rate_slow = 1000
    man.draw_trajectory(q1); man.draw_trajectory(q2)
    man.run(50000)
    print("That's all folks!")
    return man

#--------------------------------
def dipole():

    L = 8
    man = qpys.Manager(L,'Dipole in a Infinite Wall')

    w1 = man.wall(charge=(1e6)*qpys.qe,mass = 100,
                  fix=True,size=vector(L/20.,2.*L,L/20.),
                  position=vector(0,-1.*L,0),axis=vector(0,1,0))

    u1 = vector(1.,0.,0.);
    u2 = vector(0,1.,0.)
    u3 = vector(1./sqrt(2.),1./sqrt(2.),0); 
    d1 = man.dipole(distance=4,axis=u1, position=vector(-4,0,0))
    d2 = man.dipole(distance=4,axis=u2, position=vector(0,0,0))
    d3 = man.dipole(distance=4,axis=u3, position=vector(4,0,0))

    man.draw_efield(vector(-2,2,0))
    man.draw_efield(vector(2,2,0))
    ee = man.efield(vector(0,0,0))
    print('Efield ',ee, ' (N/C)')
    input('What are the forces? Enter key.')
    man.clear()
    man.draw_force(d1.qs[0]); man.draw_force(d2.qs[0]); man.draw_force(d3.qs[0])
    man.draw_force(d1.qs[1]); man.draw_force(d2.qs[1]); man.draw_force(d3.qs[1])

    input('What is the momentum of the forces? Enter key.')
    man.draw_torque(d1);man.draw_torque(d2);man.draw_torque(d3);

    input('How the dipoles will move? Enter key.')
    man.draw_wrotation(d1);man.draw_wrotation(d2);man.draw_wrotation(d3);
    man.run(nsteps=20000)
    print("That's all folks!")
    return man

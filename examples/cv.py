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

print ' Exit Ctr-D'

#----------------------------------------------
def con1():
#----------------------------------------------
    """ Example: Change potencial
    @ author: J.A. Hernando , date 18/12/2012
    """
    L = 8 # size of the grid
    man = qpys.Manager(L,'Charge Potencial')

    # create the charges
    q1 = man.charge(id="e+",position=vector(5.,0.,0),fix=True)  
    q2 = man.charge(id="e+",position=vector(-5.,0.,0),fix=True)  
    q3 = man.charge(id="e-",position=vector(0.,0.,0),fix=False)  
    #man.draw_force(q3);
    man.draw_velocity(q3)

    input('How the will move? Enter key')
    man.run(100)

    input('Enter to position the (-) charge')
    q3.vpy.pos = vector(0,3.,0)
    #q3.arrow_force.visible = False
    q3.velocity = vector(0,0,0)
    q3.arrow_velocity.visible = False
    #man.draw_force(q3);
    man.draw_velocity(q3)
    man.run(15800)

    input('Enter to position the (-) charge')
    q3.vpy.pos = vector(0,6.,0)
    #q3.arrow_force.visible = False
    q3.velocity = vector(0,0,0)
    q3.arrow_velocity.visible = False
    #man.draw_force(q3);
    man.draw_velocity(q3)
    
    input('How the will move? Enter key')
    man.run(15800)

    print("That's all folks!")
    return man
  

#----------------------------------------------
def con2():
#----------------------------------------------
    """ Example: Change potencial
    @ author: J.A. Hernando , date 18/12/2012
    """
    L = 8 # size of the grid
    man = qpys.Manager(L,'Charge Potencial')

    # create the charges
    q1 = man.charge(id="e+",position=vector(0.,0.,0),fix=True)  
    q2 = man.charge(id="e+",position=vector(4.,0.,0),fix=False)  
    man.fscale = 0.5*man.fscale
    man.dt= 0.5*man.dt
    #man.draw_force(q2); 
    man.draw_velocity(q2)

    input('How the will move? Enter key')
    man.run(10000)

    input('Enter to position the (-) charge')
    q2.vpy.pos =vector(2,0,0)
    q2.velocity=vector(0,0,0)
    #man.draw_force(q2); 
    man.draw_velocity(q2)
    input('How the will move? Enter key')
    man.run(10000)

    print("That's all folks!")
    return man
      

#--------------------------------------------------
def vol1():
#--------------------------------------------------
    """ Example: Change potencial
    @ author: J.A. Hernando , date 18/12/2012
    """

    # create the grid

    L = 8 # size of the grid
    man = qpys.Manager(L,'Charge Potencial')

    # create the charges
    q1 = man.charge(id="e+",position=vector(0.,0.,0),fix=True)  

    input('Enter key to show the potential')
    man.draw_vfield(factor=1.,unit=0.25)

    input('Enter to position for energy of (-) charge')
    man.draw_vfield(factor=-1.,unit=0.25)

    print("That's all folks!")
    return man

#--------------------------------------------------
def vol2():
#--------------------------------------------------
    """ Example: Coulomb charge
    @ author: J.A. Hernando , date 18/12/2012
    """

    # create the grid
    L = 8 # size of the grid
    man = qpys.Manager(L,'Chare Potencial')

    # create the charges
    q1 = man.charge(id="e-",position=vector(0.,0.,0),fix=True)  
    man.draw_vfield(factor=1.,unit=0.5)
    man.draw_vfield(factor=-1.,unit=0.5)

    print("That's all folks!")
    return man


#--------------------------------------------------
def vol3():
#--------------------------------------------------
    """ Example: Coulomb charge
    @ author: J.A. Hernando , date 18/12/2012
    """

    # create the grid
    L = 8 # size of the grid
    man = qpys.Manager(L,'Two Charges Potencial')

    # create the charges
    q1 = man.charge(id="e+",position=vector(4.,0.,0),fix=True)  
    q2 = man.charge(id="e+",position=vector(-4.,0.,0),fix=True) 
 
    input('Enter key to see energy of a (+) change')
    man.draw_vfield(factor=1.,unit=0.25)
    #man.draw_vfield(factor=-1.,unit=0.25)

    print("That's all folks!")
    return man

#--------------------------------------------------
def vol3bis():
#--------------------------------------------------
    """ Example: Coulomb charge
    @ author: J.A. Hernando , date 18/12/2012
    """

    # create the grid
    L = 8 # size of the grid
    man = qpys.Manager(L,'Two Charges Potencial')

    # create the charges
    q1 = man.charge(id="e+",position=vector(4.,0.,0),fix=True)  
    q2 = man.charge(id="e+",position=vector(-4.,0.,0),fix=True) 
 
    input('Enter key to see energy of a (-) change')
    #man.draw_vfield(factor=1.,unit=0.25)
    man.draw_vfield(factor=-1.,unit=0.25)

    print("That's all folks!")
    return man


#--------------------------------------------------
def vol4():
#--------------------------------------------------
    """ Example: Coulomb charge
    @ author: J.A. Hernando , date 18/12/2012
    """

    # create the grid
    L = 8 # size of the grid
    man = qpys.Manager(L,'Dipole Potencial')

    # create the charges
    q1 = man.charge(id="e+",position=vector(4.,0.,0),fix=True)  
    q2 = man.charge(id="e-",position=vector(-4.,0.,0),fix=True)  
    # q3 = man.charge(id="e+",position=vector(3,0,0),fix = False)
    # q3 = man.charge(id="e-",position=vector(0,6,0),fix = False)
    input('Enter key to see potential')
    man.draw_vfield(factor=1.,unit=0.25)
    #man.draw_vfield(factor=-1.,unit=0.5)

    print("That's all folks!")
    return man

#--------------------------------------------------
def vol5():
#--------------------------------------------------
    """ Example: Coulomb charge
    @ author: J.A. Hernando , date 18/12/2012
    """

    # create the grid
    L = 8 # size of the grid
    man = qpys.Manager(L,'Trap Potencial')

    # create the charges
    d = 5
    q1 = man.charge(id="e+",position=vector(0,d,0),fix=True)
    q2 = man.charge(id="e+",position=vector(0,-1.*d,0),fix=True)
    q3 = man.charge(id="e-",position=vector(d,0,0),fix = True)
    #q4 = man.charge(id="e-",position=vector(-1.*d,0,0),fix = False)
    man.draw_vfield(factor=1.,unit=0.25)
    man.draw_vfield(factor=-1.,unit=0.25)

    print("That's all folks!")
    return man

#--------------------------------------------------
def vol6():
#--------------------------------------------------
    """ Example: Coulomb charge
    @ author: J.A. Hernando , date 18/12/2012
    """

    # create the grid
    L = 8 # size of the grid
    man = qpys.Manager(L,'Wall Potencial')

    w1 = man.wall(charge=1e-6, mass = 100, fix=True, 
                  position=vector(-8,0,0),axis=vector(1,0,0),
                  size=vector(0.4,16.,16.))


    input('Enter key to show the potencial')
    #q3 = man.charge(id="e-",position=vector(-1.*d,0,0),fix = False)
    man.draw_vfield(factor=1.,unit=0.25)
    # man.draw_vfield(factor=-1.,unit=0.25)

    print("That's all folks!")
    return man

#--------------------------------------------------
def ene1():
#--------------------------------------------------
    """ Example: Coulomb charge
    @ author: J.A. Hernando , date 18/12/2012
    """

    # create the grid
    L = 8 # size of the grid
    man = qpys.Manager(L,'Trap')

    # create the charges
    q1 = man.charge('e+',position=vector(-6,0,0,),fix=True)
    q2 = man.charge('e+',position=vector(+6,0,0,),fix=True)
    q3 = man.charge('e+',position=vector(+4,0,0,),fix=False)
    
    input('Enter key to see energy on charge q3')
    man.draw_vfield_on(q3,unit=0.25)
    man.draw_energy(q3); man.draw_velocity(q3); man.draw_force(q3)

    print('Charge in the middle is free. The other two are fix')
    input('How the will move? Enter key')

    man.run(50000)

    print("That's all folks!")
    return man

#--------------------------------------------------
def ene2():
#--------------------------------------------------
    """ Example: Coulomb charge
    @ author: J.A. Hernando , date 18/12/2012
    """

    # create the grid
    L = 8 # size of the grid
    man = qpys.Manager(L,'Trap')

    # create the charges
    q1 = man.charge('e+',position=vector(-6,0,0,),fix=True)
    q2 = man.charge('e+',position=vector(+6,0,0,),fix=True)
    q3 = man.charge('e-',position=vector(0,6,0,),fix=False)
    
    input('Enter key to see energy on charge q3 (-)')
    man.dt = 0.7*man.dt
    man.draw_vfield_on(q3,unit=0.5)
    man.draw_energy(q3); man.draw_velocity(q3); man.draw_force(q3)

    input('How it will move the negative charge? Entery key')

    man.run(50000)

    print("That's all folks!")
    return man

#--------------------------------------------------
def ene3():
#--------------------------------------------------
    """ Example: Coulomb charge
    @ author: J.A. Hernando , date 18/12/2012
    """

    # create the grid
    L = 8 # size of the grid
    man = qpys.Manager(L)

    # create the charges
    d = 5
    q1 = man.charge('e+',position=vector(0,d,0),fix=True)
    q2 = man.charge('e+',position=vector(0,-1.*d,0),fix=True)
    q3 = man.charge('e-',position=vector(d,0,0),fix=True)
    q4 = man.charge('e-',position=vector(-1.*d,0,0),fix=False)

    input('Enter key to see V')
    man.fscale = 0.5*man.fscale
    man.draw_vfield_on(q4,unit=0.25)
    man.draw_energy(q4); man.draw_velocity(q4); man.draw_force(q4)

    print('Charge on the left is free')
    input('How it will move? Entery key')

    man.run(50000)

    print("That's all folks!")
    return man


#--------------------------------------------------
def ene4():
#--------------------------------------------------
    """ Example: Coulomb charge
    @ author: J.A. Hernando , date 18/12/2012
    """

    # create the grid
    L = 8 # size of the grid
    man = qpys.Manager(L,'Wall')
    
    q1 = man.charge('e-',position=vector(-6,2,0,),fix=False)
    q2 = man.charge('e-',position=vector(-4,0,0,),fix=False)
    q3 = man.charge('e-',position=vector(-2,-2,0,),fix=False)
    
    w1 = man.wall(charge=1e-6, mass = 100, fix=True, 
                  position=vector(8,0,0),axis=vector(-1,0,0),
                  size=vector(0.4,16.,16.))

    input('Enter to see energy on charges')
    man.draw_vfield(factor=-1.,unit=0.5)

    man.draw_force(q1); man.draw_force(q2); man.draw_force(q3);
    man.draw_velocity(q1); man.draw_velocity(q2); man.draw_velocity(q3);
    man.draw_energy(q1); man.draw_energy(q2); man.draw_energy(q3)

    def bounce():
        for q in [q1,q2,q3]:
            if (q.position().x >= w1.position().x):
                q.velocity = -1.*q.velocity
            
    input('How the charges will move? Enter key')
    man.dt = man.dt*0.40
    vals = man.run(nsteps=50000,fun=bounce,rate_fun=1)
    print("That's all folks!")
    return man

#--------------------------------------------------
def ene5():
#--------------------------------------------------
    """ Example: Coulomb charge
    @ author: J.A. Hernando , date 18/12/2012
    """

    # create the grid
    L = 8 # size of the grid
    man = qpys.Manager(L,'Wall')
    
    q1 = man.charge('e-',position=vector(6,2,0,),
                    velocity=vector(-0.6e8,0,0),fix=False)
    q2 = man.charge('e-',position=vector(5,0,0,),
                    velocity=vector(-0.5e8,0,0),fix=False)
    q3 = man.charge('e-',position=vector(4,-2,0,),
                    velocity=vector(-0.25e8,0,0),fix=False)
    
    w1 = man.wall(charge=1e-6, mass = 100, fix=True, 
                  position=vector(8,0,0),axis=vector(-1,0,0),
                  size=vector(0.4,16.,16.))

    man.draw_vfield_on(q1,unit=0.5); 
    man.draw_energy(q1); man.draw_energy(q2); man.draw_energy(q3)
    man.draw_velocity(q1); man.draw_velocity(q2); man.draw_velocity(q3) 
    #man.draw_force(q1); man.draw_force(q2); man.draw_force(q3); 

    def bounce():
        for q in [q1,q2,q3]:
            if (q.position().x >= w1.position().x):
                q.velocity = -1.*q.velocity
            
    input('How the charges will move? Enter key')
    man.dt = man.dt*0.40
    #man.rate_slow = 5000
    vals = man.run(nsteps=50000,fun=bounce,rate_fun=1)
    print("That's all folks!")
    return man

#-------------------------------------------------
def ene6():
    """ Sets few charges in a constant electic field.
    Set the value of the charge for both to reach at the same time to the plane
    """
    L = 8
    man = qpys.Manager(L,'Kepler')
    
    f = qpys.me/qpys.mp # ratio of masses e-/proton
    q1 = man.charge('e-',position=vector(0,-7,0,),velocity=vector(-4.,0,0),fix=False)
    q2 = man.charge('p',position=vector(0,0,0,),velocity=f*vector(4.,0,0),fix=False)

    man.fscale = 0.5*man.fscale
    man.draw_trajectory(q1)
    man.draw_force(q1)
    man.draw_velocity(q1)
    man.draw_vfield_on(q1,unit=0.25)
    man.draw_energy(q1)

    input('How they will move? Enter key')
    man.rate_slow = 1000    
    tup = man.run(55000)

    print("That's all folks!")
    return man

#--------------------------------------------------------
def ene7():
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

    #ars1 = list(map(lambda q: man.draw_force(q),qs))
    ars2 = list(map(lambda q: man.draw_energy(q),qs))
    ars3 = list(map(lambda q: man.draw_velocity(q),qs))
    
    comment = "The negative charges will move. How? Enter key"
    input(comment)

    #man.rate_slow = 1000 # sloww the movement...
    tup = man.run(nsteps=30000)

    print("That's all folks!")
    return man




#--------------------------------------------------
def nl():
#--------------------------------------------------
    """ Example: Change potencial
    @ author: J.A. Hernando , date 18/12/2012
    """

    # create the grid

    L = 8 # size of the grid
    man = qpys.Manager(L,'Charge Potencial')

    d = 3
    # create the charges
    q1 = man.charge(id="e+",position=vector(-d,d,0),fix=True)  
    input('Enter key to show the potential')
    man.uscale = 2233540961.
    ss = man.draw_vfield(factor=1.,unit=0.25)

    input('Enter to add new charge')
    for si in ss: si.visible = False
    q2 = man.charge(id="e+",position=vector(d,d,0),fix=True)  
    input('Enter key to show the potential')
    man.uscale = 2233540961.
    ss = man.draw_vfield(factor=1.,unit=0.25)

    input('Enter to add new charge')
    for si in ss: si.visible = False
    q3 = man.charge(id="e+",position=vector(d,-d,0),fix=True)  
    input('Enter key to show the potential')
    man.uscale = 2233540961.
    ss = man.draw_vfield(factor=1.,unit=0.25)

    input('Enter to add new charge')
    for si in ss: si.visible = False
    q4 = man.charge(id="e+",position=vector(-d,-d,0),fix=True)  
    man.uscale = 2233540961.
    input('Enter key to show the potential')
    ss = man.draw_vfield(factor=1.,unit=0.25)


    print("That's all folks!")
    return man

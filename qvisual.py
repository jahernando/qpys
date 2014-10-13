from visual import scene
#from visual.graph import *
from functools import *
#from math import *

DEBUG = False

#-----------------

def clear(obs):
    """ the objs will not be visible anymore """
    for q in obs: q.visible = False

#------- mouse click

def nothing() : return True
def anothing(arg) : return True

def click(operation = nothing):
    """ Wait for a left click, apply operation, return value
    """
    look = True
    while look:
       if scene.mouse.events:
            evt = scene.mouse.getevent()
            if (evt.click == 'left'):
                val = operation()
                return val

def click_objects(obs, operation = anothing):
    """ Wait for click, if obj is selected in the list of objs 
    apply operation on it, return value
    """
    look = True
    while look:
       if scene.mouse.events:
            evt = scene.mouse.getevent()
            if (evt.click == 'left'):
                ok = False
                for ob in obs:
                    if (evt.pick == ob):
                        look = False
                        val = operation(ob)
                        return val

def click_object(ob, operation = nothing):
    """ Wait for click, if obj is selected apply operation on it, return value
    """
    look = True
    while look:
       if scene.mouse.events:
            evt = scene.mouse.getevent()
            if (evt.click == 'left'):
                ok = False
                if (evt.pick == ob):
                    look = False
                    val = operation()
                    return val

def click_position(operation = anothing):
    """ wait for click, apply operation on the click-position, return value
    """
    look = True
    while look:
       if scene.mouse.events:
            evt = scene.mouse.getevent()
            if (evt.click == 'left'):
                pos = evt.pickpos
                #print('click at ',pos)
                val = operation(pos)
                look = False
                return val

def drag_object(vpy,operation=nothing):
    """ Move the object (vpy) while drag, apply operation
    """
    obj = None
    look = True
    while look:
        if scene.mouse.events:
            evt = scene.mouse.getevent()
            if (evt.click == 'left'):
                #print('Click, returning')
                look = False
                return
            elif ((evt.drag) and (evt.pick == vpy)):
                dpos = evt.pickpos
                obj = evt.pick
            elif (evt.drop):
                obj = None
        if (obj):
            npos = scene.mouse.project(normal=(0,0,1))
            obj.pos += npos - dpos
            dpos = npos
            # print("position ",obj.pos)
            val = operation()
            #q.set_position(obj.pos)
            #step(qs)
                
def drag_rotate(vpy,operation=nothing):
    obj = None
    look = True
    while look:
        if scene.mouse.events:
            evt = scene.mouse.getevent()
            if (evt.click == 'left'):
                #print('Click, returning')
                look = False
                return
            elif ((evt.drag) and (evt.pick == vpy)):
                dpos = evt.pickpos
                obj = evt.pick
            elif (evt.drop):
                obj = None
        if (obj):
            npos = scene.mouse.project(normal=(0,0,1))
            vpy.axis = npos
            #print("axis ",npos)
            operation()
            #q.set_position(obj.pos)
            #step(qs)

def pause():
    """ Pause for a left click of the mouse of for a keyboard input
    """
    while True:
        rate(50)
        if scene.mouse.events:
            m = scene.mouse.getevent()
            if m.click == 'left': return
        elif scene.kb.keys:
            k = scene.kb.getkey()
            return


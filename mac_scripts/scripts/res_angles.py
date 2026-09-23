import numpy as np

def min180To180(val):
    while val < -np.pi:
        val += 2*np.pi
    while val > np.pi:
        val -= 2*np.pi
    return (val*180/np.pi)
vmin180To180=np.vectorize(min180To180)

def res21(w1,w2,l1,l2):
    theta1 = vmin180To180( 2*l2 - l1 - w1 )
    theta2 = vmin180To180( 2*l2 - l1 - w2 )
    return theta1, theta2

def res32(w1,w2,l1,l2):
    theta1 = vmin180To180( 3*l2 - 2*l1 - w1 )
    theta2 = vmin180To180( 3*l2 - 2*l1 - w2 )
    return theta1, theta2

def res43(w1,w2,l1,l2):
    theta1 = vmin180To180( 4*l2 - 3*l1 - w1 )
    theta2 = vmin180To180( 4*l2 - 3*l1 - w2 )
    return theta1, theta2

def res53(w1,w2,l1,l2):
    theta1 = vmin180To180( 5*l2 - 3*l1 - 2*w1 )
    theta2 = vmin180To180( 5*l2 - 3*l1 - 2*w2 )
    theta3 = vmin180To180( 5*l2 - 3*l1 - w1 - w2 )
    return theta1, theta2, theta3

def res85(w1,w2,l1,l2):
    theta1 = vmin180To180( 8*l2 - 5*l1 - 3*w1 )
    theta2 = vmin180To180( 8*l2 - 5*l1 - 3*w2 )
    theta3 = vmin180To180( 8*l2 - 5*l1 - 2*w1 - w2 )
    theta4 = vmin180To180( 8*l2 - 5*l1 - w1 - 2*w2 )
    return theta1, theta2, theta3, theta4

def res121(l1,l2,l3):
    return vmin180To180( l1 - 2*l2 + l3)

def res132(l1,l2,l3):
    return vmin180To180( l1 - 3*l2 + 2*l3)
        
def res253(l1,l2,l3):
    return vmin180To180( 2*l1 - 5*l2 + 3*l3)


# coding: utf-8

# In[23]:


import numpy as np
import pandas as pd
import math
import serial
import time
ser = serial.Serial('/dev/ttyACM0', 9600)

def send_data (t1servo,t2servo,t3servo,t4servo) :
    data=str(t1servo)+'a'+str(t2servo)+'a'+str(t3servo)+'a'+str(t4servo)
    print(data)
    b = bytes(data, 'utf-8')
    ser.write(b)

def anglecalculator(x,y,z,t1,t2):
    x0=(a2+a3)*np.sin(t2)*np.sin(t1)
    y0=(a2+a3)*np.sin(t2)*np.cos(t1)
    z0=(a2+a3)*np.cos(t2)
    zbar=z-a1
    r1=np.array([x0 , y0 , z0])
    r2=np.array([x-x0 , y-y0 , zbar-z0])
    t4=np.arccos(np.dot(r1,r2)/((a2+a3)*a4))  
    if np.isnan(t4)==False and t4>=(5/180)*np.pi and t4<=(77/180)*np.pi:
        a=(a4*np.cos(t4))/(a2+a3)
        xc=(a+1)*x0
        yc=(a+1)*y0
        zc=(a+1)*z0
        r3=np.array([xc-x0 , yc-y0 , zc-z0])
        xe=((a2+a3)*np.sin(t2) + a4*np.cos(t4))*np.sin(t1)
        ye=((a2+a3)*np.sin(t2) + a4*np.cos(t4))*np.cos(t1)
        ze=(a2+a3)*np.cos(t2) - a4*np.sin(t4)
        r4=np.array([x-xc , y-yc , zbar-zc])
        r5=np.array([xe-xc , ye-yc , ze-zc])
        t3=np.arccos((np.dot(r4,r5))/(a4*np.sin(t4)*a4*np.sin(t4))) 
        if np.isnan(t3)==False :#and t3>=(-90/180)*np.pi and t3<=(90/180)*np.pi:
            rt=np.array([xe-x0 , ye-y0 , ze-z0])
            # n = r3 X rt
            n=[(yc*ze - yc*z0 - y0*ze - ye*zc + ye*z0 + y0*zc),
              (xe*zc - xe*z0 - x0*zc - xc*ze + xc*z0 + x0*ze),
              (xc*ye - xc*y0 - x0*ye - xe*yc - xe*y0 - x0*yc)]
            if np.dot(n,r4)<=0:
                t3=-t3
        else:
            t3=10000
    else:
        t3=10000
        t4=10000
    return ((t3/np.pi)*180,(t4/np.pi)*180)   # function returns t3,t4 in degrees

def forwardkinematics (t1test,t2test,t3test,t4test):
    # convert to radians
    t1test=(t1test/180.00)*np.pi
    t2test=(t2test/180.00)*np.pi
    t3test=(t3test/180.00)*np.pi
    t4test=(t4test/180.00)*np.pi
    #parameter table
    pt=[[(90/180)*np.pi+t1test,(90/180)*np.pi,0,a1],
        [(90/180)*np.pi+t2test,(90/180)*np.pi,0,0],
        [t3test,(-90/180)*np.pi,0,a2+a3],
        [t4test,0,a4,0]]

    i=0
    h0_1=[[np.cos(pt[i][0]) , -np.sin(pt[i][0])*np.cos(pt[i][1]) , np.sin(pt[i][0])*np.sin(pt[i][1]) , pt[i][2]*np.cos(pt[i][0])],
          [np.sin(pt[i][0]) , np.cos(pt[i][0])*np.cos(pt[i][1]) , -np.cos(pt[i][0])*np.sin(pt[i][1]) , pt[i][2]*np.sin(pt[i][0])],
          [0 , np.sin(pt[i][1]) , np.cos(pt[i][1]) , pt[i][3]],
          [0 , 0 , 0 , 1]]
    i=1
    h1_2=[[np.cos(pt[i][0]) , -np.sin(pt[i][0])*np.cos(pt[i][1]) , np.sin(pt[i][0])*np.sin(pt[i][1]) , pt[i][2]*np.cos(pt[i][0])],
          [np.sin(pt[i][0]) , np.cos(pt[i][0])*np.cos(pt[i][1]) , -np.cos(pt[i][0])*np.sin(pt[i][1]) , pt[i][2]*np.sin(pt[i][0])],
          [0 , np.sin(pt[i][1]) , np.cos(pt[i][1]) , pt[i][3]],
          [0 , 0 , 0 , 1]]
    i=2
    h2_3=[[np.cos(pt[i][0]) , -np.sin(pt[i][0])*np.cos(pt[i][1]) , np.sin(pt[i][0])*np.sin(pt[i][1]) , pt[i][2]*np.cos(pt[i][0])],
          [np.sin(pt[i][0]) , np.cos(pt[i][0])*np.cos(pt[i][1]) , -np.cos(pt[i][0])*np.sin(pt[i][1]) , pt[i][2]*np.sin(pt[i][0])],
          [0 , np.sin(pt[i][1]) , np.cos(pt[i][1]) , pt[i][3]],
          [0 , 0 , 0 , 1]]
    i=3
    h3_4=[[np.cos(pt[i][0]) , -np.sin(pt[i][0])*np.cos(pt[i][1]) , np.sin(pt[i][0])*np.sin(pt[i][1]) , pt[i][2]*np.cos(pt[i][0])],
          [np.sin(pt[i][0]) , np.cos(pt[i][0])*np.cos(pt[i][1]) , -np.cos(pt[i][0])*np.sin(pt[i][1]) , pt[i][2]*np.sin(pt[i][0])],
          [0 , np.sin(pt[i][1]) , np.cos(pt[i][1]) , pt[i][3]],
          [0 , 0 , 0 , 1]]
    h0_2=np.dot(h0_1,h1_2)
    h2_4=np.dot(h2_3,h3_4)
    h0_4=np.dot(h0_2,h2_4)
    x,y,z=h0_4[:3,3]
    return (x,y,z)

   


# In[24]:


#inputs   actual rov coordinates
#x= -15.60#543#6
#y= 887.25#773#13
#z= 416.27#-162#20

#link length
#x= 7  #6
#y= 10  #13
#z= 25  #20

#actual rov link lenghts
#a1=0 
#a2=525#529.7
#a3=90
#a4=483#411.41
#4DOF test model link lenghts
a1=12.8
a2=9.8
a3=4.4
a4=5


def inversekinematics (x,y,z):
    #define the required lists
    t1t2=[]
    angles=[]
    ang1=[]
    d0_4=[]
    angcoord=[]
    error=[]

    #create a list containing (t1,t2)
    for i in np.arange(-18,36,.5):
        for j in np.arange(29,101,1):
            t1t2.append((i,j))

    #create new list containing (t1,t2,t3,t4)
    for i in range(len(t1t2)):  
        angles.append(t1t2[i] + anglecalculator(x , y , z , (t1t2[i][0]/180)*np.pi , (t1t2[i][1]/180)*np.pi))
    
    #remove all null angles from list   
    for i in range(len(angles)):
        if (angles[i][2]!=572957.7951308233 and angles[i][3]!=572957.7951308233):   #this value is radian of 10000
            ang1.append(angles[i])
    #print(ang1) 

    #check forward kinematics
    for i in range(len(ang1)):
        t1test=-ang1[i][0]
        t2test=90 - ang1[i][1]
        t3test=ang1[i][2]
        t4test=-(90 + ang1[i][3])  
        #  actual arm mapping to t4test -(180-armangle +90)  
        # t4 mapping to actual arm 180-t4
        angcoord.append(ang1[i] + forwardkinematics(t1test,t2test,t3test,t4test))
    
    for i in range(len(angcoord)):
        error.append(math.sqrt((angcoord[i][4]-x)**2 + (angcoord[i][5]-y)**2 + (angcoord[i][6]-z)**2))

    df=pd.DataFrame(angcoord, error,columns=['t1','t2','t3','t4','x','y','z'])
    df=df.reset_index()
    df.rename(columns={'index':'error'}, inplace=True)
    #df=df.loc[(df['t3']>89) & (df['t3']<91) | (df['t3']<-89) & (df['t3']>-91) | (df['t3']>-1) & (df['t3']<1)  | (df['t3']>179) & (df['t3']<181) | (df['t3']<-179) & (df['t3']>-181)]
    #df=df.loc[(df['t1']>-1) & (df['t1']<1) | (df['t3']>179) & (df['t3']<181) | (df['t3']<-179) & (df['t3']>-181)]
    if not df.empty :
        dfmin=df.loc[df['error'].idxmin()]
    else :
        dfmin=df#.fillna(0)
    print(x,y,z)
    #print(dfmin)
    #print(df)
    return dfmin

dfmin=inversekinematics(7,10,25)
print(dfmin)

#mapping calculated angles to 4DOF test model angles
t1servo=-3.39*(dfmin['t1']-35)
t2servo=130-0.89*(dfmin['t2'])
t3servo=95+0.84*(dfmin['t3'])
t4servo=85+1.05*(dfmin['t4'])
send_data(t1servo,t2servo,t3servo,t4servo)


# In[6]:


import pygame
import time

axisvalue = np.array([0.0,0.0,0.0,0.0])
del_xyz= np.array([0.0,0.0,0.0])
x=7
y=10
z=25

pygame.init() 
#Loop until the user clicks the close button.
done = False
# Used to manage how fast the screen updates
clock = pygame.time.Clock()
# Initialize the joysticks
pygame.joystick.init()
#select and initialize joystick from which data is being retreived
joystick = pygame.joystick.Joystick(0)
joystick.init()

# -------- Main Program Loop -----------
while done==False: 
    # EVENT PROCESSING STEP
    for event in pygame.event.get(): # User did something
        if(event.type == pygame.JOYHATMOTION): # If user clicked close
            done=True # Flag that we are done so we exit this loop
        
        # Possible joystick actions: JOYAXISMOTION JOYBALLMOTION JOYBUTTONDOWN JOYBUTTONUP JOYHATMOTION                 
        if(event.type == pygame.JOYBUTTONDOWN):
            for i in range( 4 ):  # joystick has 4 axes 0 to 3
                axis = joystick.get_axis( i )
                axisvalue[i]= float("{:>6.3f}".format(axis))
            if(joystick.get_button(0)==True):
                del_xyz[0]=axisvalue[0]#*10
                del_xyz[1]=axisvalue[1]#*10
            if(joystick.get_button(1)==True):
                del_xyz[2]=axisvalue[1]#*10
            #pygame.event.set_blocked(pygame.JOYAXISMOTION)
            print(del_xyz)
            tempx=x
            tempy=y
            tempz=z
            x+=del_xyz[0]
            y+=del_xyz[1]
            z+=del_xyz[2]
	    del_xyz[0]=0
	    del_xyz[1]=0
	    del_xyz[2]=0
            dfmin=inversekinematics(x , y , z)
            if dfmin.empty :
                x=tempx
                y=tempy
                z=tempz
                dfmin=inversekinematics(x , y , z)
            print(dfmin)
            
            #mapping calculated angles to 4DOF test model angles
            t1servo=-3.39*(dfmin['t1']-35)
            t2servo=130-0.89*(dfmin['t2'])
            t3servo=95+0.84*(dfmin['t3'])
            t4servo=85+1.05*(dfmin['t4'])
            send_data(t1servo,t2servo,t3servo,t4servo)
            #pygame.event.set_allowed(pygame.JOYAXISMOTION)
    
    # Limit to 20 frames per second
    clock.tick(5)
    
# Close the window and quit.
# If you forget this line, the program will 'hang'
# on exit if running from IDLE.
pygame.quit ()
    



import matplotlib.animation as anim
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
import numpy as np
import os


def readGridFile(file:str):
    f = open(file,'r')
    gridsize = int(f.readlines(1)[0])
    vals = ' '.join(f.readlines())
    vals.replace('\n',' ')
    f.close()
    
    data = np.fromstring(vals,sep=' ').reshape(gridsize,gridsize,gridsize)
    return data,gridsize

def plotGrid(gridarray,gridsize):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    
    x,y,z,c = [],[],[],[]
    for index,val in np.ndenumerate(gridarray):
        if abs(val) > 1e-1: 
            xx,yy,zz = index
            x.append(xx - gridsize/2)
            y.append(yy - gridsize/2)
            z.append(zz - gridsize/2)
            c.append(val)
    ax.scatter(x,y,z,c=c,cmap='coolwarm')
    
    plt.show()
    plt.close()

fig = plt.figure()
fig.suptitle("This animation is dedicated to Rebecca Moore.\n She's not dead but goddamn",size='x-large')
gs = fig.add_gridspec(ncols=2)
ax1 = fig.add_subplot(gs[:])
# ax1 = fig.add_subplot(gs[:],projection="3d")
# ax2 = fig.add_subplot(gs[1],projection="3d")
arrowScaleFactor = 3

def update3D(frame):
    ax1.clear()
    # ax2.clear()
    # for a in [ax1,ax2]:
    #     a.set_xlim(-6,6)
    #     a.set_ylim(-6,6)
    #     a.set_zlim(-6,6)
    ex = f"outputs/ex/ex{frame}.dat"
    ey = f"outputs/ey/ey{frame}.dat"
    ez = f"outputs/ez/ez{frame}.dat"
    hx = f"outputs/hx/hx{frame}.dat"
    hy = f"outputs/hy/hy{frame}.dat"
    hz = f"outputs/hz/hz{frame}.dat"
    jx = f"outputs/jx/jx{frame}.dat"
    jy = f"outputs/jy/jy{frame}.dat"
    p = f"outputs/p/p{frame}.dat"

    data1,size1 = readGridFile(hx)
    data2,size2 = readGridFile(hy)
    data3,size3 = readGridFile(hz)
    
    x,y,z = [],[],[]
    u,v,w = [],[],[]

    cols = []
    cmap = cm.coolwarm
    norm = Normalize()
    norm.autoscale(cols)
    for index,val in np.ndenumerate(data1):
        xx,yy,zz = index
        x.append(xx - size1/2)
        y.append(yy - size1/2)
        z.append(zz - size1/2)
        if abs(val) > 1e+1: 
            u.append(val/arrowScaleFactor)
        else:
            u.append(0)
    for index,val in np.ndenumerate(data2):
        if abs(val) > 1e+1:
            v.append(val/arrowScaleFactor)
        else: v.append(0)
    for index,val in np.ndenumerate(data3):
        if abs(val) > 1e+1:
            w.append(val/arrowScaleFactor)
            cols.append(val)
        else: 
            w.append(val)
            cols.append(val)
    ax1.quiver(x,y,z,u,v,w,colors=cmap(norm(cols)),normalize=False)

    # data1,size1 = readGridFile(jx)
    # data2,size2 = readGridFile(jy)

    # x,y,z = [],[],[]
    # c = []
    # for index,val in np.ndenumerate(data1):
    #     if abs(val) > 1e-3:
    #         xx,yy,zz = index
    #         x.append(xx - size1/2)
    #         y.append(yy - size1/2)
    #         z.append(zz - size1/2)
    #         c.append(val)
    # ax1.scatter(x,y,z,c=c)

    # x,y,z = [],[],[]
    # c = []
    # for index,val in np.ndenumerate(data2):
    #     if abs(val) > 1e-2:
    #         xx,yy,zz = index
    #         x.append(xx - size2/2)
    #         y.append(yy - size2/2)
    #         z.append(zz - size2/2)
    #         c.append(val)
    # ax2.scatter(x,y,z,c=c)

# hz = f"outputs/hz/hz10.dat"
# data,size = readGridFile(hz)
# data = data.reshape(size,size,size)
# ax1.imshow(data[:,:,size//2])
# print(np.min(data),np.max(data))

def update2D(frame):
    ax1.clear()

    hz = f"outputs/hz/hz{frame}.dat"
    data,size = readGridFile(hz)
    data = data.reshape(size,size,size)
    ax1.imshow(data[:,:,size//2],vmin=-1,vmax=1,cmap='coolwarm')
    # print(min(data),max(data))

ani = anim.FuncAnimation(fig=fig,func=update2D,frames=200,interval=50)
ani.save(filename="figs/test_hfield.gif", writer="pillow")

plt.show()


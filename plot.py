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
gs = fig.add_gridspec(ncols=3)
ax1 = fig.add_subplot(gs[:])
# ax2 = fig.add_subplot(gs[1])
# ax3 = fig.add_subplot(gs[2])
# ax1 = fig.add_subplot(gs[0],projection="3d")
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

    data1,size1 = readGridFile(ex)
    data2,size2 = readGridFile(ey)
    data3,size3 = readGridFile(ez)
    
    xe,ye,ze = [],[],[]
    ue,ve,we = [],[],[]

    cols = []
    cmap = cm.coolwarm
    norm = Normalize()
    norm.autoscale(cols)
    for index,val in np.ndenumerate(data1):
        xx,yy,zz = index
        xe.append(xx - size1/2)
        ye.append(yy - size1/2)
        ze.append(zz - size1/2)
        if abs(val) > 1e-5: 
            ue.append(val/arrowScaleFactor)
        else:
            ue.append(0)
    for index,val in np.ndenumerate(data2):
        if abs(val) > 1e-5:
            ve.append(val/arrowScaleFactor)
        else: ve.append(0)
    for index,val in np.ndenumerate(data3):
        if abs(val) > 1e-5:
            we.append(val/arrowScaleFactor)
            cols.append(val)
        else: 
            we.append(0)
            cols.append(val)
    ax1.quiver(xe,ye,ze,ue,ve,we,colors=cmap(norm(cols)),normalize=True)

    ax2.clear()
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
        if abs(val) > 1e-5: 
            u.append(val/arrowScaleFactor)
        else:
            u.append(0)
    for index,val in np.ndenumerate(data2):
        if abs(val) > 1e-5:
            v.append(val/arrowScaleFactor)
        else: v.append(0)
    for index,val in np.ndenumerate(data3):
        if abs(val) > 1e-5:
            w.append(val/arrowScaleFactor)
            cols.append(val)
        else: 
            w.append(0)
            cols.append(val)
    ax2.quiver(x,y,z,u,v,w,colors=cmap(norm(cols)),normalize=True)


    for i in range(len(w)):
        vec1 = np.array([u[i],v[i],w[i]])
        vec2 = np.array([ue[i],ve[i],we[i]])
        print(np.dot(vec1,vec2))

def update2D(frame):
    ax1.clear()
    # Variable plot
    var = "hz"
    p = f"outputs/{var}/{var}{frame}.dat"
    p2 = f"outputs/hcon/hcon{frame}.dat"
    data1,size = readGridFile(p)
    data2,_ = readGridFile(p2)
    # d = np.divide(data2,data1)
    d = data1.reshape(size,size,size)
    m = 1e-2 #np.max(np.abs(data1))
    im = ax1.imshow(d[:,:,size//2],cmap="coolwarm",vmin=-m,vmax=m)

    ax1.set_title(r"$H_z$")



def update2Darrow(frame):
    ax1.clear()

    ex = f"outputs/hx/hx{frame}.dat"
    ey = f"outputs/hy/hy{frame}.dat"
    datax,size = readGridFile(ex)
    datay,_ = readGridFile(ey)
    # print("Total charge: ",np.sum(np.abs(data)))
    datax = datax.reshape(size,size,size)
    datay = datay.reshape(size,size,size)

    x,y,u,v = [],[],[],[]
    for index,val in np.ndenumerate(datax):
        xx,yy,_ = index
        x.append(xx)
        y.append(yy)
        u.append(val)
    for index,val in np.ndenumerate(datay):
        v.append(val)

    side = (size-1)//2
    ax1.quiver(x,y,u,v,cmap='coolwarm')#,extent=[-side,side,-side,side])
    
    print(np.min(datax),np.max(datax))

# update2D(450)
# update2Darrow(0)
# update3D(5)

ani = anim.FuncAnimation(fig=fig,func=update2D,frames=300,interval=50)
ani.save(filename="figs/test.gif", writer="pillow")

plt.show()


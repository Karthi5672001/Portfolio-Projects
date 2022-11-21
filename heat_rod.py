import netgen.gui
from ngsolve import Draw, Redraw

from netgen.csg import *

cy = Cylinder ( Pnt(0, 0.5, 0.5), Pnt(1, 0.5, 0.5), 0.2)

left  = Plane (Pnt(0,0,0), Vec(-1,0,0) ).bc("left")
right = Plane (Pnt(1,1,1), Vec( 1,0,0) ).bc("right")
front = Plane (Pnt(0,0,0), Vec(0,-1,0) ).bc("neumann")
back  = Plane (Pnt(1,1,1), Vec(0, 1,0) ).bc("neumann")
bot   = Plane (Pnt(0,0,0), Vec(0,0,-1) ).bc("left")
top   = Plane (Pnt(1,1,1), Vec(0,0, 1) ).bc("right")

cube = left * right * front * back * bot * top

geo = CSGeometry()
geo.Add (cy*cube)
ngmesh = geo.GenerateMesh(maxh=0.5)
mesh = Mesh(ngmesh)
fes = H1(mesh, order=1, dirichlet="left|right")

time = 0.0
dt = 0.001
tstep = 2 
kappa = 0.1

# define trial- and test-functions
u = fes.TrialFunction()
v = fes.TestFunction()

# the bilinear-form 
a = BilinearForm(fes, symmetric=True)
a += kappa * grad(u)*grad(v)*dx
a.Assemble()

m = BilinearForm(fes, symmetric=False)
m += u*v*dx
m.Assemble()

# the right hand side
f = LinearForm(fes)
f += 0 * v * dx

mstar = m.mat.CreateMatrix()
mstar.AsVector().data = m.mat.AsVector() + dt * a.mat.AsVector()

gfu = GridFunction(fes)
gfu.vec.data[:]=0.0
gfu.Set(x, definedon=mesh.Boundaries("left"))
gfu.Set(x, definedon=mesh.Boundaries("right"))
Draw(gfu,mesh,"u")


res = f.vec.CreateVector()
invmstar = mstar.Inverse(freedofs=fes.FreeDofs())

t_intermediate=0 # time counter within one block-run
while t_intermediate < tstep - 0.5 * dt:
    res.data = dt * f.vec - dt * a.mat * gfu.vec
    gfu.vec.data += invmstar * res
    t_intermediate += dt
    print("\r",time+t_intermediate,end="")
    Redraw(blocking=True)
print("")
time+=t_intermediate

# plot the solution (netgen-gui only)
Draw (gfu)
Draw (-grad(gfu), mesh, "Flux")

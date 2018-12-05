 #ES EDITABLE EN ANDROID CON LA APLICACIÓN SIMPLE TEXT EDITOR
println("============================================")
println("ATRACTOR DE ROSSLER")
println("Por: Dr. Jorge Chávez Carlos ")
println("============================================")
println("El modelo de Rossler se esta resolviendo...")
using ODE;
function f(t, r)
    # Extract the coordinates from the r vector
    (x, y, z) = r
    
    # The Rossler equations
    dx_dt = -y - z
    dy_dt = x + a*y
    dz_dt = b+ z*(x - c)
    
    # Return the derivatives as a vector
    [dx_dt; dy_dt; dz_dt]
end;
# Define time vector and interval grid
const dt = 0.001
const tf = 100.0
t = 0:dt:tf

# Initial position in space
const r0 = [0.1; 0.0; 0.0]

# Constants sigma, rho and beta
const a = 1.0/4.0
const b = 1.0
const c = 5.5;
(t, pos) = ode45(f, r0, t)
x = map(v -> v[1], pos)
y = map(v -> v[2], pos)
z = map(v -> v[3], pos);
println("La solucion se ha encontrado")
using PyPlot
plot3D(x, y, z);
fig, ax = subplots(1, 3, sharex=true, sharey=true, figsize=(16,8))

ax[1][:plot](x, y)
ax[1][:set_title]("X-Y cut")

ax[2][:plot](x, z)
ax[2][:set_title]("X-Z cut")

ax[3][:plot](y, z)
ax[3][:set_title]("Y-Z cut");

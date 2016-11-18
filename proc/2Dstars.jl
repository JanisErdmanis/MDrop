using PyPlot

phi = 0:0.1:(2*pi)
n = 3

r = 1 + 0.1*(cos(3*phi) + cos(2*phi))

x = r.*cos(phi)
y = r.*sin(phi)

plot(x,y)

### rinkis

r2 = 3

x2 = r2.*cos(phi)
y2 = 2*r2.*sin(phi)

plot(x2,y2)

### rinka perturbeesana
### perturbacijas lielumam jabuut proporcionalam cos(phi) kvadratam

### Citaadi perturbeet var, ja phi nolasa no normales virziena
### Var perturbeet vienaliciigi ar visaam iespeejamaam perturbaacijaam, kur maksimalaa frekvence ir limiteeta ar rezga izmeeru

phi3 = atan2(y2,x2)
r3 = sqrt(x2.^2 + y2.^2) + 0.1*cos(phi3*n)

x3 = r3.*cos(phi3)
y3 = r3.*sin(phi3)
plot(x3,y3)

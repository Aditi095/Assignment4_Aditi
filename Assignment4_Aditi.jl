
using LinearAlgebra, Statistics;
function Newtonfixedpointmap(f, f_prime, x_0, tolerance, maxiter)
    # setup the algorithm
    x_old = x_0
    normdiff = Inf
    iter = 1
    while normdiff > tolerance && iter <= maxiter
        x_new = x_old - (f(x_old)/ f_prime(x_old)) # use the passed in map
        normdiff = norm(x_new - x_old)
        x_old = x_new
        iter = iter + 1
    end
    return (x_old, normdiff, iter)
end

maxiter = 1000
tolerance = 1.0E-7
x_0 = 0.0 # initial condition

f(v) = (v-1)^3 
f_prime(v)= 3*(v-1)^2 

x_star, normdiff, iter = Newtonfixedpointmap(f, f_prime, x_0, tolerance, maxiter)
println("Fixed point = $x_star, and |f(x) - x| = $normdiff in $iter iterations")

f(v)= 27*v^3 - 3*v + 1
f_prime(v)= 81*v^2 - 3

x_star, normdiff, iter = Newtonfixedpointmap(f, f_prime, x_0, tolerance, maxiter)
println("Fixed point = $x_star, and |f(x) - x| = $normdiff in $iter iterations")

using ForwardDiff

# operator to get the derivative of this function using AD
D(f) = v-> ForwardDiff.derivative(f, v)

function Newtonfixedpointmap(f, x_0, tolerance, maxiter)
    # setup the algorithm
    x_old = x_0
    normdiff = Inf
    iter = 1
    f_prime = D(f)
    while normdiff > tolerance && iter <= maxiter
        x_new = x_old - (f(x_old)/ f_prime(x_old)) # use the passed in map
        normdiff = norm(x_new - x_old)
        x_old = x_new
        iter = iter + 1
    end
    return (x_old, normdiff, iter)
end
maxiter = 1000
tolerance = 1.0E-7
x_0 = 0.0 # initial condition

f(v) = (v-1)^3 
f_prime(v)= 3*(v-1)^2 

x_star, normdiff, iter = Newtonfixedpointmap(f, f_prime, x_0, tolerance, maxiter)
println("Fixed point = $x_star, and |f(x) - x| = $normdiff in $iter iterations")

f(v)= 27*v^3 - 3*v + 1
f_prime(v)= 81*v^2 - 3

x_star, normdiff, iter = Newtonfixedpointmap(f, f_prime, x_0, tolerance, maxiter)
println("Fixed point = $x_star, and |f(x) - x| = $normdiff in $iter iterations")

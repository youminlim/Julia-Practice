# Applies Runge Kutta method to the system

function rkSolver(rk, f, t0, tf, y0, n)
    # Declare arrays
    h = (tf - t0) / n
    t = range(t0, tf, n)
    y = zeros(n, length(y0))

    # Add in initial conditions
    y[1, :] = y0

    # RK order checker
    if rk == 1
        for i in range(1, n-1)
            k1 = f(t[i], y[i, :])
            y[i+1, :] = y[i, :] + (h * (k1/6))
        end
    elseif rk == 2
        for i in range(1, n-1)
            k1 = f(t[i], y[i, :])
            k2 = f(t[i] + h, y[i, :] + (h*k1))
            y[i+1, :] = y[i, :] + (h * (k1/2 + k2/2))
        end
    elseif rk == 3
        for i in range(1, n-1)
            k1 = f(t[i], y[i, :])
            k2 = f(t[i] + (0.5*h), y[i, :] + (0.5*h*k1))
            k3 = f(t[i] + h, y[i, :] + (h*(-k1+2*k2)))
            y[i+1, :] = y[i, :] + (h * (k1/6 + k2/3 + k3/6))
        end
    elseif rk == 4
        for i in range(1, n-1)
            k1 = f(t[i], y[i, :])
            k2 = f(t[i] + (0.5*h), y[i, :] + (0.5*h*k1))
            k3 = f(t[i] + (0.5*h), y[i, :] + (0.5*h*k2))
            k4 = f(t[i] + h, y[i, :] + (h*k3))
            y[i+1, :] = y[i, :] + (h * (k1/6 + k2/3 + k3/3 + k4/6))
        end
    else
        printIn("Order of Runge-Kutta method not defined")
    end

    return t, y
end

function twoBody3D(t, y)
    return [y[7], y[8], y[9], y[10], y[11], y[12], (G*m2*(y[4] - y[1])/r^3), (G*m2*(y[5] - y[2])/r^3), (G*m2*(y[6] - y[3])/r^3), (G*m2*(y[1] - y[4])/r^3), (G*m2*(y[2] - y[5])/r^3), (G*m2*(y[3] - y[6])/r^3)]
end

function eulerian(vector1, vector2)
    # Calculates the Eulerian distance between two vectors
    distance = vector2 - vector1
    if length(distance) == 1
        return distance
    else
        squares = 0
        for i in distance
            squares += i^2
        end
        return sqrt(squares)
    end
end

function plotter(t, y)
    plt = plot3d
end

# System parameters
G = 6.67 * 10 ^-11
m1 = 1000
m2 = 10
y0 = [0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0] # In the form y0 = [x01 y01 z01 x02 y02 z02 dot(x01 y01 z01 x02 y02 z02)]
r = eulerian(y0[1:3], y0[4:6])

rk = 4 # Choice between 1, 2, 3 and 4

# Time interval
t0 = 0      # Initial time
tf = 10    # Final time
n = 11

t, y = rkSolver(rk, twoBody3D, t0, tf, y0, n)

for i in range(1, length = n) 
    time = t[i]
    position = y[i]
    println("$time, $position")
end
# Functions to solve and visualise 2 Body Problem

function twoBody3D(t, y, r)
    return [y[7], y[8], y[9], y[10], y[11], y[12], (G*m2*(y[4] - y[1])/r^3), (G*m2*(y[5] - y[2])/r^3), (G*m2*(y[6] - y[3])/r^3), (G*m1*(y[1] - y[4])/r^3), (G*m1*(y[2] - y[5])/r^3), (G*m1*(y[3] - y[6])/r^3)]
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


function plot2Body(t, y)
    body1 = y[:, 1:3]
    body2 = y[:, 4:6]

    plot(body1[:, 1], body1[:, 2], body1[:, 3], title = "2 Body Problem", aspect_ratio = :equal, label = "Mass 1")
    plot!(body2[:, 1], body2[:, 2], body2[:, 3], label = "Mass 2")
end


function rkSolver(rk, f, t0, tf, y0, n)
    # Declare arrays
    h = (tf - t0) / n
    t = range(t0, tf, n)
    y = zeros(n, length(y0))

    # Add in initial conditions
    y[1, :] = y0
    r = eulerian(y[1, 1:3], y[1, 4:6])

    # RK order checker
    if rk == 1
        for i in range(1, n-1)
            println(r)
            k1 = f(t[i], y[i, :], r)
            y[i+1, :] = y[i, :] + (h * (k1/6))
            r = eulerian(y[i+1, 1:3], y[i+1, 4:6])
        end
    elseif rk == 2
        for i in range(1, n-1)
            k1 = f(t[i], y[i, :], r)
            k2 = f(t[i] + h, y[i, :] + (h*k1), r)
            y[i+1, :] = y[i, :] + (h * (k1/2 + k2/2))
            r = eulerian(y[i+1, 1:3], y[i+1, 4:6])
        end
    elseif rk == 3
        for i in range(1, n-1)
            k1 = f(t[i], y[i, :], r)
            k2 = f(t[i] + (0.5*h), y[i, :] + (0.5*h*k1), r)
            k3 = f(t[i] + h, y[i, :] + (h*(-k1+2*k2)), r)
            y[i+1, :] = y[i, :] + (h * (k1/6 + k2/3 + k3/6))
            r = eulerian(y[i+1, 1:3], y[i+1, 4:6])
        end
    elseif rk == 4
        for i in range(1, n-1)
            k1 = f(t[i], y[i, :], r)
            k2 = f(t[i] + (0.5*h), y[i, :] + (0.5*h*k1), r)
            k3 = f(t[i] + (0.5*h), y[i, :] + (0.5*h*k2), r)
            k4 = f(t[i] + h, y[i, :] + (h*k3), r)
            y[i+1, :] = y[i, :] + (h * (k1/6 + k2/3 + k3/3 + k4/6))
            r = eulerian(y[i+1, 1:3], y[i+1, 4:6])
        end
    else
        printIn("Order of Runge-Kutta method not defined")
    end

    return t, y
end


# Physics Checker
function gravity(G, m, y)
    # Calculates the gravitational force experienced by a body

    r = eulerian(y[1:3], y[4:6])
    return (G * m) / r^2
end


function accelChecker(G, m1, m2, y)
    # Checks if the acceleration forces experienced by the bodies are of correct magnitude
    
    dims = size(y)
    for i in range(1, length = dims[1])
        if gravity(G, m1, y[i, :]) != -(m1 / m2) * gravity(G, m2, y[i, :])
         println("Error!")
        end 
     end 
end
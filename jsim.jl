using GLMakie, BenchmarkTools, Printf

const LENGTH::Int64 = 100 # width of the lattice
const STRENGTH::Float64 = 1.0 # >1 -> ferromagnetic, <1 -> antiferromagnetic
const EXTERNAL::Float64 = 0.0 # external magnetic field
const ITERATIONS::Int64 = 10^7 # how long to run the model for

# define a structure for the spin system
mutable struct SpinSystem
    spins::Array{Int64, 2} 
    J::Float64 # interaction strength between lattice sites
    h::Float64 # external magnetic field
    L::Int64 # size of the lattice
end

# define a function to initialize the spin system
function init_spin_system(L::Int64, J::Float64, h::Float64)
    # initialize the system with up spins and down spins 
    spins = rand([-1, 1], L, L)
    return SpinSystem(spins, J, h, L)
end


function get_neighbors(i::Int64, j::Int64, L::Int64)
    # get the neighbors of the lattice site (i, j)
    neighbors = [(i, mod1(j + 1, L)), (i, mod1(j - 1, L)), (mod1(i + 1, L), j), (mod1(i - 1, L), j)]
    return neighbors
end

# function that computes the energy of a given spin system configuration
function compute_energy(system::SpinSystem)
    energy = 0.0
    # loop over all lattice sites
    for i in 1:system.L
        for j in 1:system.L
            for neighbor in get_neighbors(i, j, system.L)
                energy += -system.J * system.spins[i, j] * system.spins[neighbor[1], neighbor[2]]
            end
            # external field energy
            energy += -system.h * system.spins[i, j]
        end
    end 

    return energy
end

function delta_e(system::SpinSystem, i::Int64, j::Int64)
    # compute the change in energy if the spin at (i, j) is flipped
    delta = 0.0
    for neighbor in get_neighbors(i, j, system.L)
        delta += -system.J * system.spins[i,j]*(-system.spins[neighbor[1], neighbor[2]])
    end
    delta += -2*system.h * system.spins[i, j]
end

# computes the net magnetization of the system
function magnetization(system::SpinSystem)
    return sum(system.spins)
end

# computes the net magnetization of the system
function magnetization_density(system::SpinSystem)
    return magnetization(system) / (system.L^2)
end


# function that plots the spin system
function plot_spin_system(system::SpinSystem)
    # create a new figure
    fig = Figure(resolution = (800, 800))
    # create a new axis
    ax = Axis(fig[1, 1])
    # plot the spins
    heatmap!(ax, system.spins, markersize = 20, color = :black, colormap= :julia)
    display(fig)
end


# the main loop
function main()
    # create a spin system 
    system = init_spin_system(LENGTH,STRENGTH,EXTERNAL)
    # make an observable for the spins
    s_obs = Node(system.spins)
    # create the figure and axis
    fig = Figure(resolution = (800, 800))
    ax = Axis(fig[1, 1])
    # everytime the observable updates, we replot the spins 
    heatmap!(ax, s_obs, markersize = 20, color = :black, colormap= :julia)
    display(fig)

    # for loop that just updates the spins
    for iteration in 1:ITERATIONS
        # generate two random numbers between 1 and L
        i = rand(1:system.L)
        j = rand(1:system.L)

        # compute the change in energy if the spin at (i, j) is flipped
        delta = delta_e(system, i, j)

        # if the change in energy is negative, flip the spin
        if delta < 0
            system.spins[i, j] *= -1
        # if the change in energy is positive, flip the spin with probability exp(-delta_e)
        elseif rand() < exp(-delta)
            system.spins[i, j] *= -1
        end
        s_obs[] = system.spins
    end
end

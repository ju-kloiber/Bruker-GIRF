using BrukerGIRF
using GLMakie   # for plotting

# folder containing the simulation files
folder = "./resources/SQ2"

# some parameters
tDwell      = 8e-3      # dwell time for the gradients, in ms
RampTime    = 0.110     # time for the gradient ramps, in ms

# load gradients
t, G = load_gradients(folder, tDwell, RampTime)

# plot gradients
fig = Figure()
ax = Axis(fig[1, 1], xlabel="Time (ms)", ylabel="G (mT/m)")
lines!(ax, t, G[:, 1], color=:red, label="Gx")
lines!(ax, t, G[:, 2], color=:green, label="Gy")
lines!(ax, t, G[:, 3], color=:blue, label="Gz")
axislegend(ax)
xlims!(ax, (0, 50))
save(joinpath("images", "gradients_plot.png"), fig)
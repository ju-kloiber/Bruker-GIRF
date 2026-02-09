using BrukerGIRF
using GLMakie   # for plotting

# folder containing the simulation files
folder = "./resources/SQ2"

# some parameters
tDwell      = 8e-3      # dwell time for the gradients, in ms
RampTime    = 0.110     # time for the gradient ramps, in ms

# load gradients
G, t = load_gradients(folder, tDwell, RampTime)

# load GIRF from .mat file
matfile = joinpath("output", "exp7", "girf_filtered.mat")
girf_dict = matread(matfile)
girfT = girf_dict["girfT"]
t_girf = girf_dict["t_meas"]

# apply GIRF to gradients
G_corrected, t_corrected = apply_girf(G, t, girfT, t_girf, :direct)

# plot gradients
fig = Figure()
ax = Axis(fig[1, 1], xlabel="Time (ms)", ylabel="G (mT/m)")
lines!(ax, t, G[:, 1], color=:red, label="Gx")
lines!(ax, t, G[:, 2], color=:green, label="Gy")
lines!(ax, t, G[:, 3], color=:blue, label="Gz")
axislegend(ax)
xlims!(ax, (0, 50))
save(joinpath("images", "gradients.png"), fig)

fig = Figure()
ax = Axis(fig[1, 1], xlabel="Time (ms)", ylabel="G (mT/m)")
lines!(ax, t, G[:, 1], color=:black, label="Gx nominal", linestyle=:dash)
lines!(ax, t, G_corrected[:, 1], color=:red, label="Gx filtered")
axislegend(ax)
xlims!(ax, (0, 50))
save(joinpath("images", "gradients_raw_vs_corrected.png"), fig)
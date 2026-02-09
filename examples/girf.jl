using BrukerGIRF
using MAT
using FFTW
using GLMakie


### SOME PARAMETERS ###
dict_girf = Dict{String, Any}()    # initialize dictionary to hold parameters

# name of the folder in /resources
dict_girf["experiment"] = "exp7"
dict_girf["folder"] = "resources" 

# list of measurement files
dict_girf["files"] = ["result_7.jcamp", "result_14.jcamp", "result_21.jcamp", "result_29.jcamp", "result_36.jcamp", "result_43.jcamp", "result_50.jcamp", 
         "result_58.jcamp", "result_65.jcamp", "result_72.jcamp", "result_80.jcamp", "result_87.jcamp", 
         "result_94.jcamp"]

# corresponding gradient amplitudes in mT/m
dict_girf["amps"] = [48.601, 97.201, 145.802, 194.402, 243.003, 291.603, 340.204, 388.804, 437.405, 486.005, 534.606, 583.207, 631.807]

# some hardcoded parameters
dict_girf["p"] = 6075.067308     # max slew rate in mT/m/ms
dict_girf["tDwell"] = 2e-3       # dwell time of data acquisition in ms
dict_girf["nPoints"] = 5028      # number of points in each measurement (from Bruker header)
dict_girf["gamma"] = 42.577      # gyromagnetic ratio in MHz/T

# parameters for the raised cosine filter
dict_girf["BW"] = 60     # bandwidth in kHz
dict_girf["α"] = 1/4     # roll-off factor

### COMPUTATIONS ###
# load measured B and k from files
G, t_meas       = load_G(dict_girf)     # B from files seems to be processed already
k_meas, t_meas  = load_k(dict_girf)     # size: num_timepoints x 3 x num_amps

# compute measured gradient waveforms from k-space data
G_comp, t_meas = compute_G(k_meas, t_meas, dict_girf)

# construct nominal gradient waveforms
G_nom, t_nom = construct_nominal_G(t_meas, dict_girf)   # size: num_timepoints x 3 x num_amps

# compute GIRF
girfF, girfT, freqs, t = compute_girf(G_nom, G_comp, t_meas)

# save GIRF to .mat file
mkpath(joinpath("output", dict_girf["experiment"])) # make sure output directory exists
matfile = joinpath("output", dict_girf["experiment"], "girf_unfiltered.mat")
mat_dict = Dict("girfF" => girfF, "girfT" => girfT, "freqs" => freqs, "t_meas" => t_meas)
matwrite(matfile, mat_dict)

# construct and apply filter for GIRF
rc = raised_cosine(freqs, dict_girf)
girfF_filtered = girfF .* reshape(rc, :, 1, 1)
girfT_filtered = real.(fftshift(ifft(ifftshift(girfF_filtered, 1), 1), 1))

# save filtered GIRF to .mat file
matfile = joinpath("output", dict_girf["experiment"], "girf_filtered.mat")
mat_dict = Dict("girfF" => girfF_filtered, "girfT" => girfT_filtered, "freqs" => freqs, "t_meas" => t_meas)
matwrite(matfile, mat_dict)

# apply GIRF to nominal gradient to get corrected gradient
G_corr = apply_girf(G_nom, t_nom, girfT_filtered, t, :direct)

### PLOTS AND SAVING RESULTS ###
mkpath(joinpath("images", dict_girf["experiment"]))  # make sure image output directory exists

# plot nominal and measured Gx, Gy, Gz for each amplitude
# Gx
fig = Figure(size = (1200, 800), fontsize = 32)
ax1 = Axis(fig[1, 1], xlabel="Time (ms)", ylabel="Gradient Amplitude (mT/m)", title ="Input and Output Gradient Waveforms")
for i in eachindex(dict_girf["files"])
    lines!(ax1, t_meas, G_nom[:, 1, i], linestyle=:dash, color=:black, linewidth=5)
    lines!(ax1, t_meas, G_comp[:, 1, i], label="$(round(dict_girf["amps"][i], digits=1)) mT/m", linewidth=5)
end
axislegend(ax1)
xlims!(ax1, (-0.05, 0.3))
save(joinpath("images", dict_girf["experiment"], "measured_vs_nominal_x.png"), fig)

#Gy
fig = Figure()
ax2 = Axis(fig[1, 1], xlabel="Time (ms)", ylabel="Measured Gradient (mT/m)", title ="Gradient along Y")
for i in eachindex(dict_girf["files"])
    lines!(ax2, t_meas, G_comp[:, 2, i], label="$(round(dict_girf["amps"][i], digits=1)) mT/m", linewidth=3)
    lines!(ax2, t_meas, G_nom[:, 2, i], linestyle=:dash, color=:black, linewidth=2)
end
axislegend(ax2)
xlims!(ax2, (-0.05, 0.3))
save(joinpath("images", dict_girf["experiment"], "measured_vs_nominal_y.png"), fig)

#Gz
fig = Figure()
ax3 = Axis(fig[1, 1], xlabel="Time (ms)", ylabel="Measured Gradient (mT/m)", title ="Gradient along Z")
for i in eachindex(dict_girf["files"])
    lines!(ax3, t_meas, G_comp[:, 3, i], label="$(round(dict_girf["amps"][i], digits=1)) mT/m", linewidth=3)
    lines!(ax3, t_meas, G_nom[:, 3, i], linestyle=:dash, color=:black, linewidth=2)
end
axislegend(ax3)
xlims!(ax3, (-0.05, 0.3))
save(joinpath("images", dict_girf["experiment"], "measured_vs_nominal_z.png"), fig)

# plot GIRF in frequency domain
fig = Figure(size = (1000, 1200), fontsize = 32)
ax1 = Axis(fig[1, 1], xlabel="Frequency (kHz)", ylabel="GIRF magnitude", title ="GIRF in Frequency Domain")
lines!(ax1, freqs, abs.(girfF[:, 1, 1]), linewidth=5, label = "x -> x")
ylims!(ax1, (0, 1.1))
xlims!(ax1, (-50, 50))
save(joinpath("images", dict_girf["experiment"], "girf_unfiltered.png"), fig)

# plot filtered GIRF
fig = Figure()
ax1 = Axis(fig[1, 1])
lines!(ax1, freqs, abs.(girfF[:, 1, 1]), linewidth=2, label="Raw GIRF Magnitude")
lines!(ax1, freqs, abs.(girfF_filtered[:, 1, 1]), linewidth=3, label="Filtered GIRF Magnitude")
lines!(ax1, freqs, rc, linewidth=5, linestyle = :dash, color = :black, label="Raised Cosine Filter")
axislegend(ax1)
ylims!(ax1, (0, 1.2))

# plot corrected vs measured gradient for selected amplitudes
fig = Figure(size = (1200, 1000), fontsize = 32)
ax = Axis(fig[1, 1], xlabel="Time (ms)", ylabel="Gradient (mT/m)", title ="Corrected and Measured Gradient")
lines!(ax, t_meas, G_nom[:, 1, 12], linewidth=3, linestyle=:dash, color=:black, label = "Nominal")
lines!(ax, t_meas, G_comp[:, 1, 12], linewidth=7, color = :grey60, label = "Measured")
lines!(ax, t_meas, ifftshift(G_corr[:, 1, 12]), linewidth=5, color = :red, linestyle = :dash, label = "Computed")
lines!(ax, t_meas, G_nom[:, 1, 9], linewidth=3, linestyle=:dash, color=:black, label = "Nominal")
lines!(ax, t_meas, G_comp[:, 1, 9], linewidth=7, color = :grey60, label = "Measured")
lines!(ax, t_meas, ifftshift(G_corr[:, 1, 9]), linewidth=5, color = :red, linestyle = :dash, label = "Computed")
lines!(ax, t_meas, G_nom[:, 1, 4], linewidth=3, linestyle=:dash, color=:black, label = "Nominal")
lines!(ax, t_meas, G_comp[:, 1, 4], linewidth=7, color = :grey60, label = "Measured")
lines!(ax, t_meas, ifftshift(G_corr[:, 1, 4]), linewidth=5, color = :red, linestyle = :dash, label = "Computed")
axislegend(ax, unique=true)
xlims!(ax, (-0.05, 0.25))
save(joinpath("images", dict_girf["experiment"], "corrected_vs_measured_gradient.png"), fig)
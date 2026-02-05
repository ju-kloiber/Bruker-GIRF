using DelimitedFiles
using FFTW
using DSP

# this function loads the B field data from the specified files
function load_G(params::Dict{String, Any})

    # unpack some parameters
    nPoints     = params["nPoints"]
    files       = params["files"]
    experiment  = params["experiment"]
    folder      = params["folder"]
    tDwell      = params["tDwell"]

    # preallocate array to hold B field data
    B = zeros(nPoints, 3, length(files))

    # create time array
    t = collect(0:tDwell:(nPoints-1) * tDwell)

    # load B from each file
    for (i, file) in enumerate(files)

        # load raw data
        raw = readdlm(joinpath(folder, experiment, file), ' ', skipstart = 10)
        data = vec(permutedims(raw))    # convert to 1D array
        data = filter(!isempty, data)   # remove empty entries

        # extract Bx, By, Bz
        idx = findall(x -> isa(x, AbstractString) && occursin("Bx", x), data)[1]
        idx_start = idx + 3
        len = data[idx + 1]
        Bx = data[idx_start:idx_start + len - 1]

        idx = findall(x -> isa(x, AbstractString) && occursin("By", x), data)[1]
        idx_start = idx + 3
        len = data[idx + 1]
        By = data[idx_start:idx_start + len - 1]

        idx = findall(x -> isa(x, AbstractString) && occursin("Bz", x), data)[1]
        idx_start = idx + 3
        len = data[idx + 1]
        Bz = data[idx_start:idx_start + len - 1]

        # store in array
        B[:, 1, i] .= Bx
        B[:, 2, i] .= By
        B[:, 3, i] .= Bz
    end
    return B, t
end

# this function loads the k-space data from the specified files
function load_k(params::Dict{String, Any})

    # unpack some parameters
    nPoints     = params["nPoints"]
    files       = params["files"]
    experiment  = params["experiment"]
    folder      = params["folder"]
    tDwell      = params["tDwell"]

    # preallocate array to hold k-space data
    k = zeros(nPoints, 3, length(files))

    # create time array
    t = collect(0:tDwell:(nPoints-1) * tDwell)

    # load k from each file
    for (i, file) in enumerate(files)
        raw = readdlm(joinpath(folder, experiment, file), ' ', skipstart = 10)
        data = vec(permutedims(raw))
        data = filter(!isempty, data)

        idx = findall(x -> isa(x, AbstractString) && occursin("Kx", x), data)[1]
        idx_start = idx + 3
        len = data[idx + 1]
        Kx = data[idx_start:idx_start + len - 1]

        idx = findall(x -> isa(x, AbstractString) && occursin("Ky", x), data)[1]
        idx_start = idx + 3
        len = data[idx + 1]
        Ky = data[idx_start:idx_start + len - 1]

        idx = findall(x -> isa(x, AbstractString) && occursin("Kz", x), data)[1]
        idx_start = idx + 3
        len = data[idx + 1]
        Kz = data[idx_start:idx_start + len - 1]

        k[:, 1, i] .= Kx
        k[:, 2, i] .= Ky
        k[:, 3, i] .= Kz
    end

    return k, t
end

function compute_G(k_meas::Array{Float64, 3}, t_meas::Vector{Float64}, params::Dict{String, Any})

    # unpack parameters
    gamma = params["gamma"]  # in MHz/T

    # compute G from k
    G_comp = diff(k_meas, dims = 1) ./ (diff(t_meas)) ./ (gamma) * 1e3 # in mT/m

    # adjust time vector for diff operation
    t_meas = t_meas[1:end-1]

    return G_comp, t_meas
end

function construct_nominal_G(t_meas::Vector{Float64}, params::Dict{String, Any})

    # unpack parameters
    nPoints     = params["nPoints"]
    p           = params["p"] # max slew rate in mT/m/ms
    files       = params["files"]
    amps        = params["amps"]    

    # preallocate nominal gradient array
    G_nom = zeros(nPoints - 1, 3, length(files))

    # compute nominal gradient waveforms
    for i in eachindex(files)

        T = amps[i] / p # in ms
        G_nom[:, 1, i] .= triangle(t_meas, T, amps[i])  # Bx
        G_nom[:, 2, i] .= triangle(t_meas, T, amps[i])  # By
        G_nom[:, 3, i] .= triangle(t_meas, T, amps[i])  # Bz
    end

    return G_nom, t_meas
end

# this function generates a triangular gradient waveform with the specified duration and amplitude
function triangle(t::AbstractArray, pulse_dur::AbstractFloat, amplitude::AbstractFloat)

    G_nom = zeros(length(t))
    for i in eachindex(t)
        τ = abs(t[i] - pulse_dur)
        if τ <= pulse_dur 
            G_nom[i] = amplitude * (1 - τ / pulse_dur)
        else
            G_nom[i] = 0.0
        end
    end
    return G_nom
end

# this function computes the GIRF from nominal and measured gradient waveforms
function compute_girf(Gnom::AbstractArray, Gmeas::AbstractArray, t::AbstractArray)

    freq, _, _ = time2freq(t)
    inFreq  = fftshift(fft(Gnom, 1), 1)   # FFT of input gradient
    outFreq = fftshift(fft(Gmeas, 1), 1)    # FFT of measured output gradient

    # Number of samples (time points) 
    nS = length(t)

    # Number of input channels
    nIn = size(Gnom, 2)

    # Number of output channels 
    nOut = size(Gmeas, 2)

    # Number of repeated measurements/waveforms
    nWaveforms = size(Gmeas, 3)

    # Allocate array for GIRF result (frequency response matrix)
    girf_f = zeros(ComplexF64, (nS, nOut, nIn))   # 1028 x 3 x 3

    # Compute inverse of sum of squares of input spectrum over waveforms
    # Shape: (nS × nIn x 1)
    inSOSInv = 1 ./ sum(abs.(inFreq[:, :, :]) .^ 2, dims = 3)

    # Loop over input-output channel combinations
    for iOut in 1:nOut
        for iIn in 1:nIn
            # Compute numerator: sum of outFreq * conj(inFreq) over repetitions
            # Multiply with inverse sum of squares to get LS estimate of H(f)
            #                           (nS x nWaveforms)            (nS x nWaveforms)              (nS)
            girf_f[:, iOut, iIn] .= sum(outFreq[:, iOut, :] .* conj.(inFreq[:, iIn, :]), dims=2) .* inSOSInv[:, iIn, 1]
        end
    end

    # Replace any NaN values with 0 (e.g., if division by 0 occurred)
    girf_f[isnan.(girf_f)] .= 0

    girf_t = real.(fftshift(ifft(ifftshift(girf_f, 1), 1), 1)) # Inverse FFT to get GIRF in time domain

   return girf_f, girf_t, freq, t
end

# this function converts a time vector to a frequency vector
function time2freq(t::Vector{T}) where T<:AbstractFloat
    nrs = length(t)         # Number of time samples
    dt = t[2] - t[1]        # Time step (assumed uniform)
    f_max = 1 / dt          # Nyquist frequency
    df = f_max / nrs        # Frequency resolution (Hz per bin)

    # Construct frequency vector centered at zero (for use with fftshift)
    f = ((0:nrs-1)' .- floor(Int, nrs / 2)) .* df

    return vec(f), df, f_max
end


function apply_girf(G_raw::AbstractArray, t_raw::AbstractArray, girfT::AbstractArray, t_girf::AbstractArray, method=:direct)


    # define some parameters
    nS = size(G_raw, 1)         # number of samples
    nIn = size(G_raw, 2)        # number of input channels, e.g., 3 for Gx, Gy, Gz
    nWaveforms = size(G_raw, 3) # number of waveforms

    # check if interpolation needed
    Δt_raw = t_raw[2] - t_raw[1]
    Δt_girf = t_girf[2] - t_girf[1]

    # check if time resolutions match, abs tol 0.5 μs
    if !isapprox(Δt_raw, Δt_girf; atol=5e-4)
        @info "Time resolutions do not match (Δt_raw = $Δt_raw ms, Δt_girf = $Δt_girf ms). Interpolating..."

        t_girf_interp = collect(t_girf[1]:Δt_raw:t_girf[end])
        girfT_interp = zeros(size(t_girf_interp, 1), size(girfT, 2), size(girfT, 3))

        for iIn in 1:nIn
            for iOut in 1:nIn
                spl = Spline1D(t_girf, girfT[:, iOut, iIn], k=1)
                girfT_interp[:, iOut, iIn] .= spl.(t_girf_interp)
            end
        end

        t_girf = t_girf_interp
        girfT = girfT_interp
    end

    # preallocate corrected gradient array
    G_corr = zeros(ComplexF64, size(G_raw))

    if method == :direct

        # loop over waveforms and input channels to apply inverse GIRF
        for iWaveform in 1:nWaveforms
            for iIn in 1:nIn
                G_corr[:, iIn, iWaveform] .= conv(girfT[:, iIn, iIn], G_raw[:, iIn, iWaveform], algorithm=:direct)[1:nS]
            end
        end
        
    elseif method == :fft

        # FFT the nominal input
        G_raw_fft = fftshift(fft(G_raw, 1), 1)
        girfF = fftshift(fft(girfT, 1), 1)

        # loop over waveforms and input channels to apply inverse GIRF
        for iWaveform in 1:nWaveforms
            for iIn in 1:nIn
                G_corr_f = G_raw_fft[:, iIn, iWaveform] .* girfF[:, iIn, iIn]
                G_corr[:, iIn, iWaveform] .= ifft(ifftshift(G_corr_f), 1)
            end
        end    
    else
        error("Unknown correction method: $method")
    end

    return real.(G_corr)
end

# this function constructs a raised cosine filter
function raised_cosine(f::Vector{S}, params::Dict{String,Any}) where {S<:Real}
    
    # unpack parameters
    BW = params["BW"]      # bandwidth in kHz
    α  = params["α"]       # roll-off factor

    # preallocate filter array
    filter = zeros(size(f))

    # Define the transition band limits
    f1 = (1 - α) * BW  # end of flat passband
    f2 = (1 + α) * BW  # end of transition band

    for (i, fi) in enumerate(f)
        afi = abs(fi)
        if afi <= f1
            # Flat passband region
            filter[i] = 1
        elseif f1 < afi <= f2
            # Raised cosine transition region
            filter[i] = 0.5 * (1 + cos(π * (afi - f1) / (2 * α * BW)))
        else
            # Stopband
            filter[i] = 0
        end
    end

    return filter
end
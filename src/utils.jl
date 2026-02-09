using DelimitedFiles
using FFTW
using DSP

"""
    load_G(params::Dict{String, Any})

Loads the B field data from the specified files in the parameters dictionary.

# Arguments
- `params::Dict{String, Any}`: A dictionary containing the following keys:
    - `"nPoints"`: Number of points in each measurement (from Bruker header)
    - `"files"`: List of measurement files
    - `"experiment"`: Name of the experiment (used to construct file paths)
    - `"folder"`: Name of the folder in /resources where the files are located
    - `"tDwell"`: Dwell time of data acquisition in ms

# Returns
- `B::Array{Float64, 3}`: A 3D array containing the B field data with dimensions (nPoints, 3, num_files)
- `t::Vector{Float64}`: A vector containing the time points corresponding to the B field data
"""
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
"""
    load_k(params::Dict{String, Any})

Loads the k-space data from the specified files in the parameters dictionary.

# Arguments
- `params::Dict{String, Any}`: A dictionary containing the following keys:
    - `"nPoints"`: Number of points in each measurement (from Bruker header)
    - `"files"`: List of measurement files
    - `"experiment"`: Name of the experiment (used to construct file paths)
    - `"folder"`: Name of the folder in /resources where the files are located
    - `"tDwell"`: Dwell time of data acquisition in ms

# Returns
- `k::Array{Float64, 3}`: A 3D array containing the k-space data with dimensions (nPoints, 3, num_files)
- `t::Vector{Float64}`: A vector containing the time points corresponding to the k-space data
"""
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

"""
    compute_G(k_meas::Array{Float64, 3}, t_meas::Vector{Float64}, params::Dict{String, Any})

Computes the gradient waveforms from the measured k-space data by differentiation.

# Arguments
- `k_meas::Array{Float64, 3}`: A 3D array containing the measured k-space data with dimensions (nPoints, 3, num_files)
- `t_meas::Vector{Float64}`: A vector containing the time points corresponding to the k-space data
- `params::Dict{String, Any}`: A dictionary containing the following key:
    - `"gamma"`: Gyromagnetic ratio in MHz/T

# Returns
- `G_comp::Array{Float64, 3}`: A 3D array containing the computed gradient waveforms with dimensions (nPoints-1, 3, num_files)
- `t_meas::Vector{Float64}`: A vector containing the time points corresponding to the computed gradient waveforms (adjusted for differentiation)
"""
function compute_G(k_meas::Array{Float64, 3}, t_meas::Vector{Float64}, params::Dict{String, Any})

    # unpack parameters
    gamma = params["gamma"]  # in MHz/T

    # compute G from k
    G_comp = diff(k_meas, dims = 1) ./ (diff(t_meas)) ./ (gamma) * 1e3 # in mT/m

    # adjust time vector for diff operation
    t_meas = t_meas[1:end-1]

    return G_comp, t_meas
end

"""
    construct_nominal_G(t_meas::Vector{Float64}, params::Dict{String, Any})

Constructs the nominal gradient waveforms based on the specified parameters. The nominal waveforms are assumed to be triangular pulses with durations determined by the gradient amplitudes and the maximum slew rate.

# Arguments
- `t_meas::Vector{Float64}`: A vector containing the time points corresponding to the measured data
- `params::Dict{String, Any}`: A dictionary containing the following keys:
    - `"nPoints"`: Number of points in each measurement (from Bruker header)
    - `"p"`: Maximum slew rate in mT/m/ms
    - `"files"`: List of measurement files
    - `"amps"`: Corresponding gradient amplitudes in mT/m

# Returns
- `G_nom::Array{Float64, 3}`: A 3D array containing the nominal gradient waveforms with dimensions (nPoints-1, 3, num_files)
- `t_nom::Vector{Float64}`: A vector containing the time points corresponding to the nominal gradient waveforms
"""
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

"""
    triangle(t::AbstractArray, pulse_dur::AbstractFloat, amplitude::AbstractFloat)

Generates a triangular gradient waveform.

# Arguments
- `t::AbstractArray`: Time vector
- `pulse_dur::AbstractFloat`: Duration of the triangular pulse (time to peak)
- `amplitude::AbstractFloat`: Peak amplitude of the triangular pulse

# Returns
- `G_nom::AbstractArray`: The generated triangular gradient waveform
"""
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

"""
    compute_girf(Gnom::AbstractArray, Gmeas::AbstractArray, t::AbstractArray)

Computes the Gradient Impulse Response Function (GIRF) in both frequency and time domains using least squares estimation.

# Arguments
- `Gnom::AbstractArray`: A 3D array containing the nominal gradient waveforms with dimensions (nPoints, 3, num_files)
- `Gmeas::AbstractArray`: A 3D array containing the measured gradient waveforms with dimensions (nPoints, 3, num_files)
- `t::AbstractArray`: A vector containing the time points corresponding to the gradient waveforms

# Returns
- `girf_f::AbstractArray`: A 3D array containing the GIRF in the frequency domain with dimensions (nPoints, 3, 3)
- `girf_t::AbstractArray`: A 3D array containing the GIRF in the time domain with dimensions (nPoints, 3, 3)
- `freq::Vector{Float64}`: A vector containing the frequency points corresponding to the GIRF
- `t::Vector{Float64}`: A vector containing the time points corresponding to the GIRF
"""
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

"""
    time2freq(t::Vector{T}) where T<:AbstractFloat

Converts a time vector to its corresponding frequency vector using FFT conventions.

# Arguments
- `t::Vector{T}`: A vector of time points

# Returns
- `f::Vector{T}`: A vector of frequency points corresponding to the input time vector
- `df::T`: Frequency resolution (Hz per bin)
- `f_max::T`: Maximum frequency (Nyquist frequency)
"""
function time2freq(t::Vector{T}) where T<:AbstractFloat
    nrs = length(t)         # Number of time samples
    dt = t[2] - t[1]        # Time step (assumed uniform)
    f_max = 1 / dt          # Nyquist frequency
    df = f_max / nrs        # Frequency resolution (Hz per bin)

    # Construct frequency vector centered at zero (for use with fftshift)
    f = ((0:nrs-1)' .- floor(Int, nrs / 2)) .* df

    return vec(f), df, f_max
end


"""
    apply_girf(G_raw::AbstractArray, t_raw::AbstractArray, girfT::AbstractArray, t_girf::AbstractArray, method=:direct)

Applies the Gradient Impulse Response Function (GIRF) to the raw gradient waveforms to obtain corrected gradients. The application can be done using either direct convolution or FFT-based multiplication.

# Arguments
- `G_raw::AbstractArray`: A 3D array containing the raw gradient waveforms
- `t_raw::AbstractArray`: A vector containing the time points corresponding to the raw gradient waveforms
- `girfT::AbstractArray`: A 3D array containing the GIRF in
- `t_girf::AbstractArray`: A vector containing the time points corresponding to the GIRF
- `method::Symbol`: Method for applying the GIRF, either `:direct` for direct convolution or `:fft` for FFT-based multiplication (default is `:direct`)

# Returns
- `G_corr::AbstractArray`: A 3D array containing the corrected gradient waveforms
- `t_corr::AbstractArray`: A vector containing the time points corresponding to the corrected gradient waveforms
"""
function apply_girf(G_raw::AbstractArray, t_raw::AbstractArray, girfT::AbstractArray, t_girf::AbstractArray, method=:direct)


    # define some parameters
    nS = size(G_raw, 1)         # number of samples
    nIn = size(G_raw, 2)        # number of input channels, e.g., 3 for Gx, Gy, Gz
    nWaveforms = size(G_raw, 3) # number of waveforms

    # check if interpolation needed
    Δt_raw = t_raw[2] - t_raw[1]
    Δt_girf = t_girf[2] - t_girf[1]

    # check if fft is allowed
    if method == :fft && Δt_raw != Δt_girf
        @warn "FFT-based GIRF application requires matching time resolutions. Switching to direct convolution."
        method = :direct
    end

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

        shift = floor(Int, size(girfT, 1) / 2) # amount to shift for convolution to be centered

        # loop over waveforms and input channels to apply inverse GIRF
        for iWaveform in 1:nWaveforms
            for iIn in 1:nIn
                conv_result = conv(girfT[:, iIn, iIn], G_raw[:, iIn, iWaveform], algorithm=:direct)
                G_corr[:, iIn, iWaveform] .= Δt_raw / Δt_girf *conv_result[shift:shift+nS-1]#[shift:shift+nS-1]
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

    return real.(G_corr), t_raw
end

"""
    raised_cosine(f::Vector{S}, params::Dict{String,Any}) where {S<:Real}

Constructs a raised cosine filter based on the specified bandwidth and roll-off factor.

# Arguments
- `f::Vector{S}`: A vector of frequency points
- `params::Dict{String,Any}`: A dictionary containing the following keys:
    - `"BW"`: Bandwidth in kHz
    - `"α"`: Roll-off factor

# Returns
- `filter::Vector{S}`: A vector containing the raised cosine filter values at the specified frequency points
"""
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
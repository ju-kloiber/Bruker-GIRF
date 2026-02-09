using DelimitedFiles
using DSP
using Dierckx

"""
    load_gradient_events(folder::String, params::Dict{String,Any})

Parses the pulseprogram file in the specified folder to extract gradient events.

# Arguments
- `folder::String`: Path to the folder containing the pulseprogram file.
- `params::Dict{String,Any}`: Dictionary of parameters used for evaluating delays.

# Returns
- `events_parsed::Vector{Tuple{Float64, String}}`: A vector of tuples where each tuple contains a delay in milliseconds and the corresponding command string.
"""
function load_gradient_events(folder::String, params::Dict{String,Any})

    # empty array for delay - command pairs
    events = Vector{Tuple{String,String}}()

    # extract defintions from the pulseprogram
    delay_defs = Dict{String, String}()

    # flag indicating main sequence
    in_main_seq = false

    # (nested) exec blocks will be skipped
    exec_depth = 0

    # if statements are evaluated
    if_depth = 0

    # determines whether if statement will be skipped
    skip_if  = false    

    # path to pulseprogram file
    file = joinpath(folder, "pulseprogram")

    # loop over each line of the pulseprogram
    for line in eachline(file)

        # remove leading or trailing whitespaces
        line = strip(line)

        # skip subroutine instantiations in the pulseprogram
        startswith(line, "subr ") && continue

        # skip EXEC blocks (can be nested)
        if occursin("EXEC_begin", line)
            exec_depth += 1
            continue
        end

        if occursin("EXEC_end", line)
            exec_depth -= 1
            continue
        end

        # 0 if same number of starts and ends
        exec_depth > 0 && continue

        # skip empty lines and pure comments
        isempty(line) && continue
        startswith(line, ";") && continue
        startswith(line, "#") && continue

        # skip single-line if statements with goto, TODO: implement properly
        if is_single_line_if(line)
            continue
        end

        # capture delay definitions, requires structure: define delay name
        if startswith(line, "define delay ")
            name = split(line)[3]
            delay_defs[name] = ""      # placeholder
            continue
        end

        # capture quoted delay definitions, requires structure "delay_name = expression"
        if startswith(line, "\"") && endswith(line, "\"")

            # remove quotes
            inner = strip(line[2:end-1])  

            # match name and expr
            m = match(r"^(\w+)\s*=\s*(.+)$", inner)
            if m !== nothing
                name = m.captures[1]
                expr = strip(m.captures[2])

                # add to definitions dict if key exists (assumes chronological order)
                if haskey(delay_defs, name)
                    delay_defs[name] = expr
                else
                    @warn "Delay definition for uninitialized delay name: $name"
                end
            end

            continue
        end
        
        # detect start of main sequence (keyword: "start,")
        if !in_main_seq && match(r"^start\s*,", lowercase(line)) !== nothing
            in_main_seq = true
        end

        # skip lines that are not in main sequence
        !in_main_seq && continue

        #######################################################
        # here begins the parsing of main sequence lines ######
        #######################################################

        # handle if statements
        if startswith(line, "if")
            cond = eval_if_condition(line, params)
            if_depth += 1
            skip_if = !cond
            continue
        end

        # IF block braces
        if line == "{"
            continue
        end

        if line == "}"
            if if_depth > 0
                if_depth -= 1
                if if_depth == 0
                    skip_if = false
                end
                continue
            end
        end

        skip_if && continue

        # Remove trailing comments
        line = split(line, ";", limit=2)[1] |> strip

        # Loop start: <label>, <delay> <command...>
        m = match(r"^(\w+)\s*,\s*(\S+)\s+(.*)$", line)
        if m !== nothing
            label   = m.captures[1]
            delay   = m.captures[2]     # ← keep raw
            command = m.captures[3]

            push!(events, ("0u", "loop_start:$label"))
            push!(events, (delay, command))
            continue
        end

        # Loop start without command
        m = match(r"^(\w+)\s*,\s*(\S+)$", line)
        if m !== nothing
            label = m.captures[1]
            delay = m.captures[2]

            push!(events, ("0u", "loop_start:$label"))
            push!(events, (delay, ""))
            continue
        end

        # Loop end
        m_loop_end = match(r"^lo\s+to\s+(\w+)\s+times\s+(\S+)", line)
        if m_loop_end !== nothing
            label = m_loop_end.captures[1]
            count = m_loop_end.captures[2]
            push!(events, ("0u", "loop_end:$label:$count"))
            continue
        end

        # Gradient events and gatepulses
        m_grad = match(r"^(\S+)\s+(grad_ramp\{[^}]+\}|grad_off|gatepulse)", line)
        if m_grad !== nothing
            delay   = m_grad.captures[1]
            command = m_grad.captures[2]
            push!(events, (delay, command))
            continue
        end
        
        # ADC
        if occursin("aqq", line)
            delay = "AcqTime"
            command = "ADC"
            push!(events, (delay, command))
            continue
        end

        # RF pulse events
        m_rf = match(r"^\((p\d+):[^)]*\):\S+", line)
        if m_rf !== nothing
            pulse = m_rf.captures[1]
            push!(events, (pulse, "pulse"))
            continue
        end

        # commands with spatial phase increments
        m = match(r"^(\S+)\s+([A-Za-z_]\w*\.(inc|dec))$", line)
        if m !== nothing
            delay   = m.captures[1]
            command = m.captures[2]
            push!(events, (delay, command))
            continue
        end

        # general delay + command, eg <delay> <reload B0> or <delay> <eoscnp>
        m_cmd = match(r"^(?!\()\s*(\S+)\s+([A-Za-z_]\w*(?:\s+[A-Za-z_]\w*)?)$", line)
        if m_cmd !== nothing
            delay   = m_cmd.captures[1]
            command = m_cmd.captures[2]

            push!(events, (delay, command))
            continue
        end

        # single delay, e.g., "8u" or "d10"
        m_cmd = match(r"^\s*((?:\d+(?:\.\d+)?u)|(?:d\d+))\s*$", line)
        if m_cmd !== nothing
            delay   = m_cmd.captures[1]
            command = ""
            push!(events, (delay, command))
            continue
        end
    end

    #####################################################
    # parse delays into numeric values ##################
    #####################################################
    
    events_parsed = Vector{Tuple{Float64, String}}(undef, length(events))

    for i in eachindex(events)
        delay, cmd = events[i]

        if haskey(delay_defs, delay)
            delay = delay_defs[delay]
        end

        delay = replace(delay, r"\bde\b" => "params[\"DE\"]")   # de → pre-scan delay
        # 
        delay = replace(delay, r"\bAcqTime\b" => "params[\"PVM_AcquisitionTime\"]")  # AcqTime → params[PVM_AcquisitionTime]

        # dN → d[N+1]
        delay = replace(delay, r"\bd(\d+)\b" => s"d[\1+1] * 1e3")  # s to ms
        
        # pN → p[N+1]
        delay = replace(delay, r"\bp(\d+)\b" => s"p[\1+1] * 1e-3")   # μs to ms

        # µs → ms
        delay = replace(delay, r"(\d+(?:\.\d+)?)u" => s"\1e-3")

        # ms → ms (optional)
        delay = replace(delay, r"(\d+(?:\.\d+)?)m\b" => s"\1")

        # evaluate the delay expression
        val = eval(Meta.parse(delay))

        # Store the parsed delay and command
        events_parsed[i] = (val, cmd)
    end

    return events_parsed
end

"""
    load_list(folder::String, label::String)

Loads a list of values from the acqp file in the specified folder.

# Arguments
- `folder::String`: Path to the folder containing the acqp file.
- `label::String`: The label of the list to load.

# Returns
- `list::Vector{Float32}`: A vector of Float32 values corresponding to the specified label.
"""
function load_list(folder::String, label::String)

    file = joinpath(folder, "acqp")
    raw = readdlm(file)
    data = vec(permutedims(raw))
    data = filter(!isempty, data)

    idx = findall(x -> isa(x, AbstractString) && occursin("\$$(label)=", x), data)[1]
    N = data[idx + 1]

    list = Float32.(data[idx + 3:idx+N+2])

    return list
end

"""
    load_gradient_file(filename::String)

Loads gradient labels and scaling factors from a gradient file.

# Arguments
- `filename::String`: Path to the gradient file.

# Returns
- `labels::Vector{String}`: A vector of gradient labels.
- `scales::Vector{Float64}`: A vector of scaling factors corresponding to the labels.
"""
function load_gradient_file(filename::String)
    labels = String[]
    scales = Float64[]

    for line in eachline(filename)
        line = strip(line)
        isempty(line) && continue

        # Match lines like: (0) + $g1(-1.000)
        m = match(r"\$g(\d+)\((-?[0-9.]+)\)", line)
        if m !== nothing
            push!(labels, "g" * m.captures[1])
            push!(scales, parse(Float64, m.captures[2]))
        else
            # Lines like (-0.000) → no gradient
            push!(labels, "0")
            push!(scales, 0.0)
        end
    end

    return labels, scales
end

"""
    load_parameter(folder::String, label::String, T=Float64)

Loads a single parameter from the acqp file in the specified folder.

# Arguments
- `folder::String`: Path to the folder containing the acqp file.
- `label::String`: The label of the parameter to load.
- `T`: The type to which the parameter should be converted (default: Float64).

# Returns
- `N::T`: The value of the specified parameter converted to type T.
"""
function load_parameter(folder::String, label::String, T=Float64)

    file = joinpath(folder, "acqp")
    raw = readdlm(file)
    data = vec(permutedims(raw))
    data = filter(!isempty, data)

    idx = findall(x -> isa(x, AbstractString) && occursin("\$$(label)=", x), data)[1]

    m = match(r"=(\d+)", data[idx])
    N = parse(T, m.captures[1])

    return N
end


"""
    eval_gradient_component(s::AbstractString, g, acq_spatial_phase_1, acq_spatial_phase_2)

Evaluates a gradient component expression string by substituting gradient labels and spatial phase increments.

# Arguments
- `s::AbstractString`: The gradient component expression string.
- `g`: The gradient amplitude array.
- `acq_spatial_phase_1`: The spatial phase increment array for the first dimension.
- `acq_spatial_phase_2`: The spatial phase increment array for the second dimension.

# Returns
- The evaluated gradient component value.
"""
function eval_gradient_component(s::AbstractString, g, acq_spatial_phase_1, acq_spatial_phase_2)
    
    # Replace gN with g[N+1]
    s2 = replace(s, r"g(\d+)" => x -> begin
        n = parse(Int, match(r"\d+", x).match)
        "g[$(n+1)]"
    end)
    
    # Replace ACQ_spatial_phase_1 / 2 with their numeric variable names
    s2 = replace(s2, "ACQ_spatial_phase_1" => ("(" * string(acq_spatial_phase_1) * ")"))
    s2 = replace(s2, "ACQ_spatial_phase_2" => ("(" * string(acq_spatial_phase_2) * ")"))
    
    # Evaluate the arithmetic expression
    return eval(Meta.parse(s2))
end


"""
    eval_if_condition(line::AbstractString, params::Dict{String,Any})

Evaluates a simple if condition from the pulseprogram file.

# Arguments
- `line::AbstractString`: The line containing the if condition.
- `params::Dict{String,Any}`: Dictionary of parameters used for evaluating the condition.

# Returns
- `Bool`: The result of the evaluated condition.
"""
function eval_if_condition(line::AbstractString, params::Dict{String,Any})
    # Match: if (PARAM OP VALUE)
    m = match(r"if\s*\(\s*(\w+)\s*(==|!=|<=|>=|<|>)\s*([0-9]+(?:\.[0-9]+)?)\s*\)", line)
    m === nothing && return false

    param = m.captures[1]
    op    = m.captures[2]
    value   = m.captures[3]

    if haskey(params, param)
        param = "params[\"" * param * "\"]"
    end

    if occursin("On", value)
        value = "true"
    elseif occursin("Off", value)
        value = "false"
    end

    expr_str = string(param, " ", op, " ", value)
    return eval(Meta.parse(expr_str))
end



function linearize_loops_old(events_parsed, params)

    linear = Vector{Tuple{Float64,String}}()
    loop_stack = Vector{Tuple{Int,String,Int}}()  # (start_idx, label, repeat_count)


    for (delay, cmd) in events_parsed

        if startswith(cmd, "loop_start:")
            # Push loop marker onto stack
            label = split(cmd, ":")[2]
            push!(linear, (delay, cmd))
            start_idx = length(linear)
            push!(loop_stack, (start_idx, label, 0))  # repeat_count filled at loop_end

        elseif startswith(cmd, "loop_end:")
            isempty(loop_stack) && error("loop_end without loop_start")
            parts = split(cmd, ":")
            label = parts[2]
            count_expr = parts[3]

            @assert !isempty(loop_stack) "Loop end without loop start: $label"

            start_idx, start_label, _ = pop!(loop_stack)

            @assert start_label == label "Loop label mismatch: start = $start_label, end   = $label"
            
            # Evaluate count using counters
            # Replace params in expression
            expr_str = count_expr
            for (k,v) in params
                expr_str = replace(expr_str, k => string(v))
            end
            count = eval(Meta.parse(expr_str))  # integer


            # Copy the block (excluding loop_start)
            block = linear[start_idx+1:end]
            for _ in 2:count
                append!(linear, block)
            end


            # Add the loop_end itself
            push!(linear, (delay, cmd))

        else
            push!(linear, (delay, cmd))
        end
    end

    return linear
end

"""
    linearize_loops(events_parsed, params)

Expands loops in the parsed gradient events into a linear sequence of events.

# Arguments
- `events_parsed::Vector{Tuple{Float64, String}}`: A vector of tuples where each tuple contains a delay in milliseconds and the corresponding command string.
- `params::Dict{String,Any}`: Dictionary of parameters used for evaluating loop counts.

# Returns
- `linear::Vector{Tuple{Float64, String}}`: A linearized vector of tuples with loops expanded.
"""
function linearize_loops(events_parsed, params)

    # initialize linear event list and loop stack
    global linear = Vector{Tuple{Float64,String}}()

    # loop over parsed events
    for n in eachindex(events_parsed)
        (delay, cmd) = events_parsed[n]
        if startswith(cmd, "loop_start:")  
            push!(linear, (delay, cmd))  # keep loop_start in linear list for now
        elseif startswith(cmd, "loop_end:")

            m = length(linear)

            parts = split(cmd, ":")
            label = parts[2]
            count_expr = parts[3]
            if occursin("/", count_expr)
                # Handle division in count expression
                count_expr = split(count_expr, "/")
                count_expr = Int(parse(Float64, string(params[count_expr[1]])) / parse(Float64, string(params[count_expr[2]]))) - 1
            elseif occursin("*", count_expr)
                # Handle division in count expression
                count_expr = split(count_expr, "*")
                count_expr = Int(parse(Float64, string(params[count_expr[1]])) * parse(Float64, string(params[count_expr[2]]))) - 1
            else
                count_expr = Int(parse(Float64, string(params[count_expr])))- 1
            end
            

            loop_idx = 1
            while loop_idx < m
                if startswith(linear[m - loop_idx][2], "loop_start:$label")
                    linear = vcat(linear, repeat(linear[m - loop_idx + 1:end], count_expr))
                    push!(linear, (delay, cmd))  # add the loop_end itself
                    break
                else
                    loop_idx += 1
                end
            end
        else
            push!(linear, (delay, cmd))
        end
    end

    return linear
end

"""
    load_pvm_params(folder::String)

Loads PVM parameters from the method file in the specified folder.

# Arguments
- `folder::String`: Path to the folder containing the method file.

# Returns
- `params::Dict{String,Any}`: A dictionary containing the loaded PVM parameters.
"""
function load_pvm_params(folder::String)

    file = joinpath(folder, "method")

    params = Dict{String,Any}()

    current_array_name   = nothing
    current_array_values = String[]

    for rawline in eachline(file)
        line = strip(rawline)

        isempty(line) && continue
        startswith(line, "\$\$") && continue

        # -------------------------------------------------
        # If we are currently collecting an array
        # -------------------------------------------------
        if current_array_name !== nothing
            if startswith(line, "##\$")
                # finalize array
                params[current_array_name] = parse_pvm_array(current_array_values)

                # reset state
                current_array_name   = nothing
                current_array_values = String[]
                # fall through to process this new ##$ line
            else
                append!(current_array_values, split(line))
                continue
            end
        end

        # -------------------------------------------------
        # Parameter line
        # -------------------------------------------------
        startswith(line, "##\$") || continue

        s = strip(line[4:end])   # remove ##$
        m = match(r"^(\w+)\s*=\s*(.+)$", s)
        m === nothing && continue

        name  = m.captures[1]
        value = strip(m.captures[2])

        # -------------------------------------------------
        # Array parameter
        # -------------------------------------------------
        if startswith(value, "(")
            current_array_name   = name
            current_array_values = String[]
            continue
        end

        # -------------------------------------------------
        # Scalar parameter
        # -------------------------------------------------
        params[name] = parse_pvm_scalar(value)
    end

    # finalize array if file ended unexpectedly
    if current_array_name !== nothing
        params[current_array_name] = parse_pvm_array(current_array_values)
    end

    return params
end


"""
    parse_pvm_scalar(value::AbstractString)

Parses a scalar PVM parameter value from a string.

# Arguments
- `value::AbstractString`: The string representation of the scalar value.

# Returns
- The parsed value, which can be a Bool, Int, Float64, or String.
"""
function parse_pvm_scalar(value::AbstractString)
    value == "On"  && return true
    value == "Off" && return false

    occursin(r"^\d+$", value) && return parse(Int, value)
    occursin(r"^\d+(\.\d+)?$", value) && return parse(Float64, value)

    return value
end

"""
    parse_pvm_array(tokens::AbstractArray{<:AbstractString})

Parses an array of PVM parameter values from a list of string tokens.

# Arguments
- `tokens::AbstractArray{<:AbstractString}`: An array of string tokens representing the values.

# Returns
- `out::Vector{Any}`: A vector containing the parsed values, which can be of mixed types (Int, Float64, String).
"""
function parse_pvm_array(tokens::AbstractArray{<:AbstractString})
    out = Any[]

    for t in tokens
        clean = replace(t, r"[<>]" => "")

        if occursin(r"^[+-]?\d+$", clean)
            push!(out, parse(Int, clean))
        elseif occursin(r"^[+-]?\d+(\.\d+)?$", clean)
            push!(out, parse(Float64, clean))
        else
            push!(out, clean)
        end
    end

    return out
end

is_single_line_if(line::AbstractString) = occursin(r"^\s*if\b", line) && occursin("goto", line) && occursin(";", line)

is_block_if(line::AbstractString) = occursin(r"^\s*if\s*\(", line) && !occursin("goto", line)


"""
    compute_gradients(linear::AbstractArray, g::AbstractArray, tDwell::AbstractFloat, RampTime::AbstractFloat, ACQ_spatial_phase_1::AbstractArray, ACQ_spatial_phase_2::AbstractArray)

Computes gradient waveforms from a linearized sequence of gradient events.

# Arguments
- `linear::AbstractArray`: A linearized array of gradient events (tuples of delay and command).
- `g::AbstractArray`: An array of gradient amplitudes.
- `tDwell::AbstractFloat`: The dwell time in milliseconds.
- `RampTime::AbstractFloat`: The ramp time in milliseconds.
- `ACQ_spatial_phase_1::AbstractArray`: The spatial phase increment array for the first dimension.
- `ACQ_spatial_phase_2::AbstractArray`: The spatial phase increment array for the second dimension. 

# Returns
- `t::Vector{Float64}`: A vector of time points.
- `Gx::Vector{Float64}`: A vector of gradient amplitudes in the x-direction.
- `Gy::Vector{Float64}`: A vector of gradient amplitudes in the y-direction.
- `Gz::Vector{Float64}`: A vector of gradient amplitudes in the z-direction.
"""
function compute_gradients(linear::AbstractArray, g::AbstractArray, tDwell::AbstractFloat, RampTime::AbstractFloat, ACQ_spatial_phase_1::AbstractArray, ACQ_spatial_phase_2::AbstractArray)

    # define arrays for time and gradients
    t  = Float64[0.0]
    Gx = Float64[0.0]
    Gy = Float64[0.0]
    Gz = Float64[0.0]

    # current time
    tcur = 0.0

    # indices for spatial phase increments
    spatial_phase_idx1 = 1
    spatial_phase_idx2 = 1

    # copies for storing spatial phase indices
    copy_spatial_phase_idx1 = 1
    copy_spatial_phase_idx2 = 1

    # loop over all linearized events
    for i in eachindex(linear)

        # unzip the event
        delay, cmd = linear[i]

        # Current gradient amplitude 
        Gx_curr = i == 1 ? 0.0 : Gx[end]
        Gy_curr = i == 1 ? 0.0 : Gy[end]
        Gz_curr = i == 1 ? 0.0 : Gz[end]

        # number of time steps in this block
        n = round(Int, delay / tDwell)

        if cmd == "pulse" || occursin("ADC", cmd) || cmd == ""
            Gx_tgt = Gx_curr
            Gy_tgt = Gy_curr
            Gz_tgt = Gz_curr
        elseif occursin("grad_off", cmd)
            Gx_tgt = 0.0
            Gy_tgt = 0.0
            Gz_tgt = 0.0
        elseif occursin(r"\bACQ_spatial_phase_1\.inc\b", cmd)
            spatial_phase_idx1 += 1
            Gx_tgt = Gx_curr
            Gy_tgt = Gy_curr
            Gz_tgt = Gz_curr
        elseif occursin(r"\bACQ_spatial_phase_2\.inc\b", cmd)
            spatial_phase_idx2 += 1
            Gx_tgt = Gx_curr
            Gy_tgt = Gy_curr
            Gz_tgt = Gz_curr
        elseif occursin(r"\bACQ_spatial_phase_1\.dec\b", cmd)
            spatial_phase_idx1 -= 1
            Gx_tgt = Gx_curr
            Gy_tgt = Gy_curr
            Gz_tgt = Gz_curr
        elseif occursin(r"\bACQ_spatial_phase_2\.dec\b", cmd)
            spatial_phase_idx2 -= 1
            Gx_tgt = Gx_curr
            Gy_tgt = Gy_curr
            Gz_tgt = Gz_curr
        elseif occursin(r"\bACQ_spatial_phase_1\.store\b", cmd)
            copy_spatial_phase_idx1 = copy(spatial_phase_idx1)
            Gx_tgt = Gx_curr
            Gy_tgt = Gy_curr
            Gz_tgt = Gz_curr
        elseif occursin(r"\bACQ_spatial_phase_1\.restore\b", cmd)
            spatial_phase_idx1 = copy(copy_spatial_phase_idx1)
            Gx_tgt = Gx_curr
            Gy_tgt = Gy_curr
            Gz_tgt = Gz_curr
        elseif occursin(r"\bACQ_spatial_phase_2\.store\b", cmd)
            copy_spatial_phase_idx2 = copy(spatial_phase_idx2)
            Gx_tgt = Gx_curr
            Gy_tgt = Gy_curr
            Gz_tgt = Gz_curr
        elseif occursin(r"\bACQ_spatial_phase_2\.restore\b", cmd)
            spatial_phase_idx2 = copy(copy_spatial_phase_idx2)
            Gx_tgt = Gx_curr
            Gy_tgt = Gy_curr
            Gz_tgt = Gz_curr
        elseif occursin("grad_ramp", cmd)

            # Match everything inside the braces
            m = match(r"\{(.*)\}", cmd)

            if m !== nothing
                # Extract the matched group
                inside = m.captures[1]
                
                # Split by comma and strip whitespace
                values = strip.(split(inside, ","))
            end

            if occursin("ACQ_spatial_phase_1", cmd)
                idx_wrapped1 = ((spatial_phase_idx1 - 1) % length(ACQ_spatial_phase_1)) + 1
                idx_wrapped2 = ((spatial_phase_idx2 - 1) % length(ACQ_spatial_phase_2)) + 1
                amp1 = ACQ_spatial_phase_1[idx_wrapped1]
                amp2 = ACQ_spatial_phase_2[idx_wrapped2]
            else 
                amp1 = 1.0
                amp2 = 1.0
            end

            Gx_tgt = eval_gradient_component(values[1], g, amp1, amp2)
            Gy_tgt = eval_gradient_component(values[2], g, amp1, amp2)
            Gz_tgt = eval_gradient_component(values[3], g, amp1, amp2)

        elseif occursin("loop", cmd)
            Gx_tgt = Gx_curr
            Gy_tgt = Gy_curr
            Gz_tgt = Gz_curr
        elseif occursin(".", cmd)
            Gx_tgt = Gx_curr
            Gy_tgt = Gy_curr
            Gz_tgt = Gz_curr
        else
            #@warn "Unknown gradient command: $cmd, index $i."
            Gx_tgt = Gx_curr
            Gy_tgt = Gy_curr
            Gz_tgt = Gz_curr
        end

        is_ramp = occursin("grad_ramp", cmd) || occursin("grad_off", cmd)
        n_ramp = is_ramp ? min(round(Int, RampTime / tDwell), n) : 0

        for k in 1:n
            push!(t, tcur)

            if is_ramp && k ≤ n_ramp
                α = k / n_ramp
                push!(Gx, α * (Gx_tgt - Gx_curr) + Gx_curr)
                push!(Gy, α * (Gy_tgt - Gy_curr) + Gy_curr)
                push!(Gz, α * (Gz_tgt - Gz_curr) + Gz_curr)
            else
                push!(Gx, Gx_tgt)
                push!(Gy, Gy_tgt)
                push!(Gz, Gz_tgt)
            end
            tcur += tDwell
        end

    end

    deleteat!(t, 1)
    deleteat!(Gx, 1)
    deleteat!(Gy, 1)
    deleteat!(Gz, 1)


    if params["ACQ_dim"] == 3
        amp_factors = 0.01 .* params["PVM_EffSWh"] ./ (params["PVM_Fov"] .* 42.577)
    else
        amp_factors = 0.01 .* params["PVM_EffSWh"] ./ (params["PVM_Fov"] .* 42.577)
        push!(amp_factors, 0.01 .* params["PVM_EffSWh"] ./ (params["PVM_SliceThick"] * 42.577))
    end

    Gx .*= amp_factors[1]  # in mT/m
    Gy .*= amp_factors[2] # in mT/m
    Gz .*= amp_factors[3]  # in mT/m


    return t, [Gx Gy Gz]
end

"""
    load_gradients(folder::String, tDwell::AbstractFloat, RampTime::AbstractFloat)

Loads gradient waveforms from the specified folder.

# Arguments
- `folder::String`: Path to the folder containing the gradient files.
- `tDwell::AbstractFloat`: The dwell time in milliseconds.
- `RampTime::AbstractFloat`: The ramp time in milliseconds.

# Returns
- `t::Vector{Float64}`: A vector of time points.
- `G::Matrix{Float64}`: A matrix of gradient amplitudes (columns correspond to x, y, z).
- `params::Dict{String,Any}`: A dictionary containing the loaded PVM parameters.
"""
function load_gradients(folder::String, tDwell::AbstractFloat, RampTime::AbstractFloat)

    # load the PVM parameters from the method file
    global params = load_pvm_params(folder) # global so eval() can access
    params["NECHOES"]   = load_parameter(folder, "NECHOES", Int)  # number of echoes
    params["DE"]        = load_parameter(folder, "DE") * 1e-3 # pre-scan delay, us
    params["ACQ_dim"]   = load_parameter(folder, "ACQ_dim", Int)  # acquisition dimension
    params["NA"]        = load_parameter(folder, "NA", Int)  # number of averages
    params["NR"]        = load_parameter(folder, "NR", Int)  # number of repetitions
    params["NAE"]       = load_parameter(folder, "NAE", Int)  # number of acquired echoes
    params["NSLICES"]   = load_parameter(folder, "NSLICES", Int)  # number of slices
    params["SW_h"]      = load_parameter(folder, "SW_h")  # spectral width in Hz → ms

    # first load the delay and pulse duration arrays
    global d                        = load_list(folder, "D")
    global p                        = load_list(folder, "P")
    global g                        = load_list(folder, "ACQ_gradient_amplitude")
    ACQ_size                        = load_list(folder, "ACQ_size")
    global l                        = load_list(folder, "L")
    global ACQ_spatial_phase_1      = load_list(folder, "ACQ_spatial_phase_1")
    if params["ACQ_dim"] == 3
        global ACQ_spatial_phase_2 = load_list(folder, "ACQ_spatial_phase_2")
    else
        global ACQ_spatial_phase_2 = ones(Int(ACQ_size[2]))
    end

    params["l2"] = l[3]
    params["ACQ_size[1]"] = ACQ_size[1]
    params["ACQ_size[2]"] = ACQ_size[2]


    # load the gradient waveform file
    #=
    file_gradx = joinpath(folder, "grdprog.ax")
    labels_x, scales_x = load_gradient_file(file_gradx)

    file_grady = joinpath(folder, "grdprog.ay")
    labels_y, scales_y = load_gradient_file(file_grady) 

    file_gradz = joinpath(folder, "grdprog.az")
    labels_z, scales_z = load_gradient_file(file_gradz)
    =#

    # load the list of gradient events from the pulseprogram file
    global events = load_gradient_events(folder, params)

    

    # linearize loops
    global linear = linearize_loops(events, params)
    # return events, linear
    # compute the gradients from the linearized events
    t, G = compute_gradients(linear, g, tDwell, RampTime, ACQ_spatial_phase_1, ACQ_spatial_phase_2)

    return G, t
end

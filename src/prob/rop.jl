""
function run_rop(file, model_constructor, optimizer; kwargs...)
    return _PM.run_model(file, model_constructor, optimizer, build_rop; multinetwork=true,
        ref_extensions=[_PM.ref_add_on_off_va_bounds!, ref_add_damaged_items!], kwargs...)
end


""
function build_rop(pm::_PM.AbstractPowerModel)
    for (n, network) in _PM.nws(pm)
        variable_bus_damage_indicator(pm, nw=n)
        variable_bus_voltage_damage(pm, nw=n)

        variable_branch_damage_indicator(pm, nw=n)
        _PM.variable_branch_power(pm, nw=n)

        _PM.variable_dcline_power(pm, nw=n)

        variable_storage_damage_indicator(pm, nw=n)
        variable_storage_power_mi_damage(pm, nw=n)

        variable_gen_damage_indicator(pm, nw=n)
        variable_gen_power_damage(pm, nw=n)

        _PM.variable_load_power_factor(pm, nw=n, relax=true)
        _PM.variable_shunt_admittance_factor(pm, nw=n, relax=true)

        constraint_restoration_cardinality_ub(pm, nw=n)

        constraint_model_voltage_damage(pm, nw=n)

        for i in _PM.ids(pm, :ref_buses, nw=n)
            _PM.constraint_theta_ref(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :bus, nw=n)
            constraint_bus_damage_soft(pm, i, nw=n)
            constraint_power_balance_shed(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :gen, nw=n)
            constraint_gen_damage(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :load, nw=n)
            constraint_load_damage(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :shunt, nw=n)
            constraint_shunt_damage(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :branch, nw=n)
            constraint_branch_damage(pm, i, nw=n)
            constraint_ohms_yt_from_damage(pm, i, nw=n)
            constraint_ohms_yt_to_damage(pm, i, nw=n)

            constraint_voltage_angle_difference_damage(pm, i, nw=n)

            constraint_thermal_limit_from_damage(pm, i, nw=n)
            constraint_thermal_limit_to_damage(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :dcline, nw=n)
            _PM.constraint_dcline_power_losses(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :storage, nw=n)
            constraint_storage_damage(pm, i, nw=n)
            _PM.constraint_storage_complementarity_mi(pm, i, nw=n)
            _PM.constraint_storage_losses(pm, i, nw=n)
        end
    end


    network_ids = sort(collect(_PM.nw_ids(pm)))
    n_1 = network_ids[1]
    for i in _PM.ids(pm, :storage, nw=n_1)
        _PM.constraint_storage_state(pm, i, nw=n_1)
    end

    for n_2 in network_ids[2:end]
        for i in _PM.ids(pm, :storage, nw=n_2)
            _PM.constraint_storage_state(pm, i, n_1, n_2)
        end
        for i in _PM.ids(pm, :gen, nw=n_2)
            constraint_gen_energized(pm, i, n_1, n_2)
        end
        for i in _PM.ids(pm, :bus, nw=n_2)
            constraint_bus_energized(pm, i, n_1, n_2)
        end
        for i in _PM.ids(pm, :storage, nw=n_2)
            constraint_storage_energized(pm, i, n_1, n_2)
        end
        for i in _PM.ids(pm, :branch, nw=n_2)
            constraint_branch_energized(pm, i, n_1, n_2)
        end
        for i in _PM.ids(pm, :load, nw=n_2)
            constraint_load_increasing(pm, i, n_1, n_2)
        end
        n_1 = n_2
    end

    n_final = last(network_ids)
    constraint_restore_all_items(pm, n_final)

    objective_max_load_delivered(pm)
end


## fast ROP?
""
function run_rop_fast(file, model_constructor, optimizer; kwargs...)
    return _PM.run_model(file, model_constructor, optimizer, build_rop_fast; multinetwork=true,
        ref_extensions=[_PM.ref_add_on_off_va_bounds!, ref_add_damaged_items!], kwargs...)
end


""
function build_rop_fast(pm::_PM.AbstractPowerModel)
    for (n, network) in _PM.nws(pm)
        variable_bus_damage_indicator(pm, nw=n)
        variable_bus_voltage_damage(pm, nw=n)

        variable_branch_damage_indicator(pm, nw=n)
        _PM.variable_branch_power(pm, nw=n)

        _PM.variable_dcline_power(pm, nw=n)

        variable_storage_damage_indicator(pm, nw=n)
        variable_storage_power_mi_damage(pm, nw=n)

        variable_gen_damage_indicator(pm, nw=n)
        variable_gen_power_damage(pm, nw=n)

        _PM.variable_load_power_factor(pm, nw=n, relax=true)
        _PM.variable_shunt_admittance_factor(pm, nw=n, relax=true)

        constraint_restoration_cardinality_ub(pm, nw=n)

        constraint_model_voltage_damage(pm, nw=n)

        for i in _PM.ids(pm, :ref_buses, nw=n)
            _PM.constraint_theta_ref(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :bus, nw=n)
            constraint_bus_damage_soft(pm, i, nw=n)
            constraint_power_balance_shed(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :gen, nw=n)
            constraint_gen_damage(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :load, nw=n)
            constraint_load_damage(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :shunt, nw=n)
            constraint_shunt_damage(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :branch, nw=n)
            constraint_branch_damage(pm, i, nw=n)
            constraint_ohms_yt_from_damage(pm, i, nw=n)
            constraint_ohms_yt_to_damage(pm, i, nw=n)

            constraint_voltage_angle_difference_damage(pm, i, nw=n)

            constraint_thermal_limit_from_damage(pm, i, nw=n)
            constraint_thermal_limit_to_damage(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :dcline, nw=n)
            _PM.constraint_dcline_power_losses(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :storage, nw=n)
            constraint_storage_damage(pm, i, nw=n)
            _PM.constraint_storage_complementarity_mi(pm, i, nw=n)
            _PM.constraint_storage_losses(pm, i, nw=n)
        end
    end


    network_ids = sort(collect(_PM.nw_ids(pm)))
    n_1 = network_ids[1]
    for i in _PM.ids(pm, :storage, nw=n_1)
        _PM.constraint_storage_state(pm, i, nw=n_1)
    end

    for n_2 in network_ids[2:end]
        for i in _PM.ids(pm, :storage, nw=n_2)
            _PM.constraint_storage_state(pm, i, n_1, n_2)
        end
        for i in _PM.ids(pm, :gen, nw=n_2)
            constraint_gen_energized(pm, i, n_1, n_2)
        end
        for i in _PM.ids(pm, :bus, nw=n_2)
            constraint_bus_energized(pm, i, n_1, n_2)
        end
        for i in _PM.ids(pm, :storage, nw=n_2)
            constraint_storage_energized(pm, i, n_1, n_2)
        end
        for i in _PM.ids(pm, :branch, nw=n_2)
            constraint_branch_energized(pm, i, n_1, n_2)
        end
        for i in _PM.ids(pm, :load, nw=n_2)
            constraint_load_increasing(pm, i, n_1, n_2)
        end
        n_1 = n_2
    end

    n_final = last(network_ids)
    constraint_restore_all_items(pm, n_final)

    objective_max_load_delivered_fast(pm)
end

"add penality to unrestored components"
function objective_max_load_delivered_fast(pm::_PM.AbstractPowerModel)
    nws = _PM.nw_ids(pm)

    @assert all(!_PM.ismulticonductor(pm, n) for n in nws)

    z_demand = Dict(n => _PM.var(pm, n, :z_demand) for n in nws)
    time_elapsed = Dict(n => get(_PM.ref(pm, n), :time_elapsed, 1.0) for n in nws)
    weight = maximum([get(_PM.ref(pm,n), :weight, 1.0) for n in nws])

    z_shunt = Dict(n => _PM.var(pm, n, :z_shunt) for n in nws)
    z_gen = Dict(n => _PM.var(pm, n, :z_gen) for n in nws)
    z_storage = Dict(n => _PM.var(pm, n, :z_storage) for n in nws)
    z_bus = Dict(n => _PM.var(pm, n, :z_bus) for n in nws)
    z_branch = Dict(n => _PM.var(pm, n, :z_branch) for n in nws)

    load_weight = Dict(n =>
        Dict(i => get(load, "weight", 1.0) for (i,load) in _PM.ref(pm, n, :load))
    for n in nws)

    M = Dict()
    for n in nws
        if !isempty(load_weight[n])
            M[n] = weight*maximum(abs.(values(load_weight[n])))
        else
            M[n]= weight
        end
    end

    return JuMP.@objective(pm.model, Max,
        sum(
            time_elapsed[n]*(
                sum(load_weight[n][i]*abs(load["pd"])*z_demand[n][i] for (i,load) in _PM.ref(pm, n, :load)) +
                sum(M[n]*z_bus[n][i] for (i,bus) in _PM.ref(pm, n, :bus)) +
                sum(M[n]*z_gen[n][i] for (i,gen) in _PM.ref(pm, n, :gen)) +
                sum(M[n]*z_storage[n][i] for (i,storage) in _PM.ref(pm, n, :storage)) +
                sum(M[n]*z_shunt[n][i] for (i,shunt) in _PM.ref(pm, n, :shunt)) +
                sum(M[n]*z_branch[n][i] for (i,branch) in _PM.ref(pm, n, :branch))
            )
        for n in nws)
    )
end

"add penality to unrestored components"
function objective_max_load_delivered_fast(pm::_PM.AbstractACPModel)
    nws = _PM.nw_ids(pm)

    @assert all(!_PM.ismulticonductor(pm, n) for n in nws)

    vm_vio = Dict(n => _PM.var(pm, :vm_vio, nw=n) for n in nws)
    z_demand = Dict(n => _PM.var(pm, n, :z_demand) for n in nws)
    time_elapsed = Dict(n => get(_PM.ref(pm, n), :time_elapsed, 1.0) for n in nws)
    weight = maximum([get(_PM.ref(pm,n), :weight, 1.0) for n in nws])

    load_weight = Dict(n =>
        Dict(i => get(load, "weight", 1.0) for (i,load) in _PM.ref(pm, n, :load))
    for n in nws)

    W = Dict()
    M = Dict()
    for n in nws
        if !isempty(load_weight[n])
            W[n] = 10*maximum(abs.(values(load_weight[n])))
            M[n] = weight*maximum(abs.(values(load_weight[n])))
        else
            W[n]= 10
            M[n]= weight
        end
    end

    z_shunt = Dict(n => _PM.var(pm, n, :z_shunt) for n in nws)
    z_gen = Dict(n => _PM.var(pm, n, :z_gen) for n in nws)
    z_storage = Dict(n => _PM.var(pm, n, :z_storage) for n in nws)
    z_bus = Dict(n => _PM.var(pm, n, :z_bus) for n in nws)
    z_branch = Dict(n => _PM.var(pm, n, :z_branch) for n in nws)

    return JuMP.@objective(pm.model, Max,
        sum(
            time_elapsed[n]*(
                sum(-W[n]*vm_vio[n][i] for (i,bus) in _PM.ref(pm, n, :bus)) +
                sum(load_weight[n][i]*abs(load["pd"])*z_demand[n][i] for (i,load) in _PM.ref(pm, n, :load)) +
                sum(M[n]*z_bus[n][i] for (i,bus) in _PM.ref(pm, n, :bus)) +
                sum(M[n]*z_gen[n][i] for (i,gen) in _PM.ref(pm, n, :gen)) +
                sum(M[n]*z_storage[n][i] for (i,storage) in _PM.ref(pm, n, :storage)) +
                sum(M[n]*z_shunt[n][i] for (i,shunt) in _PM.ref(pm, n, :shunt)) +
                sum(M[n]*z_branch[n][i] for (i,branch) in _PM.ref(pm, n, :branch))
            )
        for n in nws)
    )
end

function objective_max_load_delivered_fast(pm::_PM.AbstractWRModel)
    nws = _PM.nw_ids(pm)

    @assert all(!_PM.ismulticonductor(pm, n) for n in nws)

    w_vio = Dict(n => _PM.var(pm, :w_vio, nw=n) for n in nws)
    z_demand = Dict(n => _PM.var(pm, n, :z_demand) for n in nws)
    time_elapsed = Dict(n => get(_PM.ref(pm, n), :time_elapsed, 1.0) for n in nws)
    weight = maximum([get(_PM.ref(pm,n), :weight, 1.0) for n in nws])

    load_weight = Dict(n =>
        Dict(i => get(load, "weight", 1.0) for (i,load) in _PM.ref(pm, n, :load))
    for n in nws)

    W = Dict()
    M = Dict()
    for n in nws
        if !isempty(load_weight[n])
            W[n] = 10*maximum(abs.(values(load_weight[n])))
            M[n] = weight*maximum(abs.(values(load_weight[n])))
        else
            W[n]= 10
            M[n]= weight
        end
    end

    z_shunt = Dict(n => _PM.var(pm, n, :z_shunt) for n in nws)
    z_gen = Dict(n => _PM.var(pm, n, :z_gen) for n in nws)
    z_storage = Dict(n => _PM.var(pm, n, :z_storage) for n in nws)
    z_bus = Dict(n => _PM.var(pm, n, :z_bus) for n in nws)
    z_branch = Dict(n => _PM.var(pm, n, :z_branch) for n in nws)

    return JuMP.@objective(pm.model, Max,
        sum(
            time_elapsed[n]*(
                sum(-W[n]*w_vio[n][i] for (i,bus) in _PM.ref(pm, n, :bus)) +
                sum(load_weight[n][i]*abs(load["pd"])*z_demand[n][i] for (i,load) in _PM.ref(pm, n, :load)) +
                # sum(M[n]*z_bus[n][i] for (i,bus) in _PM.ref(pm, n, :bus)) +
                sum(M[n]*z_gen[n][i] for (i,gen) in _PM.ref(pm, n, :gen)) +
                sum(M[n]*z_storage[n][i] for (i,storage) in _PM.ref(pm, n, :storage)) +
                sum(M[n]*z_shunt[n][i] for (i,shunt) in _PM.ref(pm, n, :shunt))+
                sum(M[n]*z_branch[n][i] for (i,branch) in _PM.ref(pm, n, :branch))
            )
        for n in nws)
    )
end
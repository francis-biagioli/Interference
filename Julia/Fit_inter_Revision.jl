#Set Directory
dir = "~Julia/Interference/" 
cd(dir)

#Load Packages
using DataFrames, CSV, Statistics
using MCMCChains
using Distributions
using Turing
using StatsPlots
using Plots
using LsqFit
using Roots
using LambertW
using StatsBase
using GLM
using StatsModels

#Confirm functions
plotlyjs()
gr()
lambertw(1)


# ---------------------------------------------------- Set directories, load in FoRAGE Data -------------------------------

# --- Set Output directories ---
trace_dir = "~/Figures/Fit Trace Plots Revision/"
animated_surface = "~/Figures/Animated FR Surface Revision"
sd_check_prey = "~/Figures/SD_vs_Prey_off"
sd_check_pred = "~/Figures/SD_vs_Pred_Dens"
isdir(trace_dir) || mkpath(trace_dir)
isdir(animated_surface) || mkpath(animated_surface)
isdir(sd_check_prey) || mkpath(sd_check_prey)
isdir(sd_check_pred) || mkpath(sd_check_pred)


#Load in Data simulation function
include("~/simulate_data_set_INT.jl")

# --- Read in meta and curve data ---
df_meta = CSV.read("~/FoRAGE_V5_sources_and_meta_inter.csv", DataFrame)
dimensions = df_meta.Dim
data_types = df_meta.DataOrigin
Arena_size_2D = df_meta.TwoD_Arena_size_cm2
Arena_size_3D = df_meta.ThreeD_arena_size_cm3
replenished = df_meta.PreyReplaced

df_curves = CSV.read("~/FoRAGE_db_V5_Dec_20_2024_original_curves.csv", DataFrame)
density_units_2D = df_curves.TwoDDensityUnits
density_units_3D = df_curves.ThreeDDensityUnits

# -------------------------------------------------------- Define Models to Fit -----------------------------------------------

# --- Models ---
@model function fun_res_BDE_HDE(prey_offered, prey_eaten, pred_density, starting_a, starting_h, starting_w, starting_α, starting_β, lm_intercept, AS)
    a ~ truncated(Normal(starting_a[1], 10*starting_a[1]), lower = 0)
    h ~ truncated(Normal(starting_h[1], 2*starting_h[1]), lower = 0)
    w ~ truncated(Normal(starting_w[1], 10*starting_w[1]), lower = 0)

    for i in 1:length(prey_eaten)
        predicted = a*prey_offered[i]/(1 + a*h*prey_offered[i] + w * (pred_density[i]-1/AS))

        #σ = lm_intercept + starting_α * prey_offered[i] + starting_β * pred_density[i]
        σ = max(0.1, lm_intercept + starting_α * prey_offered[i] + starting_β * pred_density[i])

        prey_eaten[i] ~ truncated(Normal(predicted, σ), 0, maximum(prey_eaten))


    end
end

# MODEL 2 - Bolker Beddington-DeAngelis LambertW function
@model function fun_res_BBDEL(prey_offered, prey_eaten, pred_density, starting_a, starting_h, starting_w, starting_α, starting_β, lm_intercept, AS)
    a ~ truncated(Normal(starting_a[1], 10*starting_a[1]), lower = 0)
    h ~ truncated(Normal(starting_h[1], 2*starting_h[1]), lower = 0)
    w ~ truncated(Normal(starting_w[1], 10*starting_w[1]), lower = 0)

    for i in 1:length(prey_eaten)
        ϕ = a/(1+w*(pred_density[i]-1/AS)) # AS is arena size
        predicted = prey_offered[i]-LambertW.lambertw(ϕ*h*prey_offered[i]*exp(-ϕ*(1-h*prey_offered[i])))/(ϕ*h)

        #σ = lm_intercept + starting_α * prey_offered[i] + starting_β * pred_density[i]
        σ = max(0.1, lm_intercept + starting_α * prey_offered[i] + starting_β * pred_density[i])

        prey_eaten[i] ~ truncated(Normal(predicted, σ),0,maximum(prey_eaten))
    
    end
end

# ------------------------------- Functions Used For Plotting --------------------------------------------------

# Beddington-DeAngelis analytic prediction (vectorized)
function predict_BDE(a,h,w,N,P,AS)
    a*N / (1 + a*h*N + w*(P - 1/AS))
end

function predict_BBDEL(a,h,w,N,P,AS)

    ϕ = a/(1 + w*(P - 1/AS))

    return N -
        LambertW.lambertw(
            ϕ*h*N*exp(-ϕ*(1-h*N))
        )/(ϕ*h)

end

# ------------------------------------------------- Create empty elements to store parameter fits -------------------------------

# --- Store fitted parameters ---
fitted_a = zeros(maximum(df_meta.Inter_ID), 1)
fitted_h = zeros(maximum(df_meta.Inter_ID), 1)
fitted_w = zeros(maximum(df_meta.Inter_ID), 1)

# --- Create Vector to store prey_scale_adjust for each study
prey_scale_adjust_vec = zeros(43)

#declair number of studies, iterations, chains, and total posterior draws and build empty posterior matricies
n_studies = 43
n_draws_per_chain = 10000
n_chains = 3
n_draws = n_draws_per_chain * n_chains

W_mat = fill(NaN, n_studies, n_draws)
A_mat = fill(NaN, n_studies, n_draws)
H_mat = fill(NaN, n_studies, n_draws)


# -------------------------------------------------------- Begin Fitting Loop ------------------------------------------------

#identify problem studies to omit from fitting
studies_to_skip = [38]

# --- Fit models across studies ---
for i = 1:43

    if i in studies_to_skip
        println("Skipping study $i")

        fitted_a[i] = NaN
        fitted_h[i] = NaN
        fitted_w[i] = NaN
        prey_scale_adjust_vec[i] = NaN

        continue
    end

# ----------------------------------------------- FoRAGE Data Standardization ------------------------------------------------
    
    println(i)
    # --- initialize empty vectors ---
        prey_offered_2 = Vector{Float64}()
        prey_offered_3 = Vector{Float64}()
        indices = Vector{Int}()
        pred_density = Vector{Float64}()

    # --- read in data from meta sheet ---
        pred_numbs = df_meta.PredPerArena[findall(df_meta.Inter_ID .== i)] # number of preds per arena
        which_ID = df_meta.FoRAGE_ID[findall(df_meta.Inter_ID .== i)] #identifyer linking Inter_ID to FoRAGE ID

    # --- Pull dimension and arena sizes from metadata ---
        dimension = dimensions[findall(df_meta.Inter_ID .== i)][1]
        arena_size_2D = Arena_size_2D[findall(df_meta.Inter_ID .== i)][1]
        arena_size_3D = Arena_size_3D[findall(df_meta.Inter_ID .== i)][1]

    # --- compute AS once per study `i` ---
    if dimension == 2
        AS = arena_size_2D
    elseif dimension == 3
        AS = arena_size_3D
    elseif dimension == 2.5
        if !ismissing(arena_size_3D)
            AS = arena_size_3D^(2.5/3)
        else
            AS = arena_size_2D^(2.5/2)
        end
    end

    # --- find all rows in curves dataset with matching FoRAGE ID ---
    for j = 1:length(which_ID) #for each set of rows with unique FoRAGE IDs 
        toappend = findall(df_curves.Column1 .== which_ID[j]) # find all columns in curves dataset with the same FoRAGE ID
        append!(indices, toappend) # concatenates all rows within the curves dataset that belong to the same Inter_ID into a single indices

        # --- Compute predator density depending on arena dimension ---
        if dimension == 2
            # Keep cm² t units
            density_list = (pred_numbs[j] ./ arena_size_2D) .* ones(length(toappend), 1) # calculates the predator density for each curves dataset row for each Inter_ID, calculating the correct predator denity for each FoRAGE ID

        elseif dimension == 3
            # Keep cm³ units
            density_list = (pred_numbs[j] ./ arena_size_3D) .* ones(length(toappend), 1) #same as above but for 3D studies

        elseif dimension == 2.5
        # approximate “2.5D” arena using geometric scaling
            #Keep cm^2.5 units
            if !ismissing(arena_size_3D)
                arena_size_2p5D = (arena_size_3D)^(2.5 / 3)
                density_list = (pred_numbs[j] ./ arena_size_2p5D) .* ones(length(toappend), 1) # same as above, but converts 2D to 2.5D
            else
                arena_size_2p5D = (arena_size_2D)^(2.5 / 2)
                density_list = (pred_numbs[j] ./ arena_size_2p5D) .* ones(length(toappend), 1) # same as above, but converts 3D to 2.5D
            end
             # Diagnostic printout for studies
            println("Study $i:")
            println("  Dim = $dimension")
            println("  Pred Number = ", pred_numbs[j])
            println("  Arena size raw (cm²) = ", arena_size_2D)
            println("  Arena size raw (cm3) = ", arena_size_3D)
            println("  Density example = ", density_list[1])
            println("  ---")
        end
        append!(pred_density, density_list) #appends values from current study i to the complete list of predator density accross all looped studies 1:i in pred_density
    end

    # --- Call in data from FoRAGE curve dataset ---
        prey_offered_2 = df_curves.TwoDPreyDensity[indices]
        prey_offered_3 = df_curves.ThreeDPreyDensity[indices]
        prey_eaten = df_curves.ForagingRate[indices]
        std_err = df_curves.StandardError[indices]
        samp_size = df_curves.SampleSize[indices]

    # --- Estimate missing std errors or sample sizes ---
    if ismissing(std_err[1])
        std_err = exp(-1.8886) .* prey_eaten .^ 0.95
    end
    if ismissing(samp_size[1])
        samp_size = 3 .* ones(length(prey_eaten), 1)
    end

    #  --- Unit conversions from m to cm ---
    # m^2 to cm^2
    units_2D = density_units_2D[indices][1]
    if !ismissing(units_2D) && units_2D == "prey per m2"
        prey_offered_2 = prey_offered_2 ./ 1e4
    #    #prey_eaten_2 = prey_eaten ./ 1e4
    end
    # m^3 to cm^3
    units_3D = density_units_3D[indices][1]
    if !ismissing(units_3D) && units_3D == "prey per m3"
        prey_offered_3 = prey_offered_3 ./ 1e6
    #    #prey_eaten_3 = prey_eaten ./ 1e6
    end

    # --- Handle 2D / 2.5D / 3D data ---
    if dimension[] == 2
        prey_offered = prey_offered_2
        #prey_eaten = prey_eaten_2
    elseif dimension[] == 3
        prey_offered = prey_offered_3
        #prey_eaten = prey_eaten_3
    elseif dimension[] == 2.5
        if ismissing(prey_offered_3[1])
            prey_offered = prey_offered_2 .^ (2.5/2)
            #prey_eaten = prey_eaten_2 .^ (2.5/2)
            pred_density = pred_density .^ (2.5/2)
        else
            prey_offered = prey_offered_3 .^ (2.5/3)
            #prey_eaten = prey_eaten_3 .^ (2.5/3)
            pred_density = pred_density .^ (2.5/3)
        end
    end

    # --- Simulate if data are means ---
    data_type = data_types[findall(df_meta.Inter_ID .== i)][1]
    if data_type == "Mean"
        prey_eaten, prey_offered, pred_density = simulate_data_set_INT(length(prey_eaten), std_err, samp_size, prey_eaten, prey_offered, pred_density)
    end

    # --- Check and rescale eaten/offered ---
    # *This block ensures that the number of prey offered does not exceed the number of prey available*
    prey_scale_test = prey_eaten[prey_offered .> 0] ./ prey_offered[prey_offered .> 0]
    if sum(prey_scale_test .> 1) > 0
        prey_scale_adjust = 2 * maximum(prey_scale_test)
    elseif mean(prey_scale_test) < 1e-5
        prey_scale_adjust = 1e-2
    else
        prey_scale_adjust = 1
    end

    prey_offered = prey_offered .* prey_scale_adjust

    # --- store prey_scale_adjust for each study in empty vector ---
    prey_scale_adjust_vec[i] = prey_scale_adjust

# ------------------------------------------- Starting Values ----------------------------------------------

    # Calculate starting values
    max_a = maximum(prey_eaten ./ prey_offered)
    starting_a = max_a #*0.5
    starting_h = (1 / maximum(prey_eaten))

    # calculate w starting value with nonlinear inverse fit
        @. int_estimate(x, p) = 1/(p[1]*x) # model
        fit = curve_fit(int_estimate, pred_density, prey_eaten, [10.0])
        starting_w = fit.param[1]

    # Choose model type
    replenished = df_meta.PreyReplaced[findall(df_meta.Inter_ID .== i)][1]

#  ------------------------------------------- Emperical Calculations for SD ----------------------------------------------

    # combine empirical prey eaten measurements
    empirical_sigma_df = DataFrame(
        prey_offered = prey_offered,
        prey_eaten = prey_eaten,
        pred_density = pred_density
    )

    empirical_sigma_plot_df = combine(
        groupby(empirical_sigma_df, [:prey_offered, :pred_density]),
        :prey_eaten => std => :sd_prey_eaten
    )

    # Remove NaN values 
    no_nan = findall(i -> !any(isnan, empirical_sigma_plot_df[i,:]), 1:nrow(empirical_sigma_plot_df))
    empirical_sigma_plot_df_no_nan = empirical_sigma_plot_df[no_nan,:]

    #Run multiple regression to calculate emperical values for α and β
    alpha_beta_fit = lm(
        @formula( sd_prey_eaten ~ prey_offered + pred_density),
    empirical_sigma_plot_df_no_nan)

    # Store values for σ calulation 
    lm_intercept = coef(alpha_beta_fit)[1]
    starting_α = coef(alpha_beta_fit)[2]
    starting_β = coef(alpha_beta_fit)[3]

# -------------------------------------------------------- Model Selection ---------------------------------------------------

    if replenished == "Y"
        model_type2 = fun_res_BDE_HDE(prey_offered, prey_eaten, pred_density, starting_a, starting_h, starting_w, starting_α, starting_β, lm_intercept, AS)
    else
        model_type2 = fun_res_BBDEL(prey_offered, prey_eaten, pred_density, starting_a, starting_h, starting_w, starting_α, starting_β, lm_intercept, AS)
    end

# -------------------------------------------------------- Model Fitting -----------------------------------------------------

    # --- Fit model ---
    fr_chain = sample(
        model_type2,
        NUTS(2000, 0.8), #num of warmup iterations and target acceptance rate
        MCMCThreads(),
        chain_type = MCMCChains.Chains,
        n_draws_per_chain, # number of sampling iterations
        initial_params = [(
        a = starting_a, # inital value for a
        h = starting_h, # inital value for h
        w = starting_w) #initial value for w
        for _ in 1:n_chains],
        n_chains # num of chains
    )

# ------------------------------------------------------ Store Posteriors ----------------------------------------------------

    #save posterior matricies
    df_chain = DataFrame(fr_chain)

    #pull posterior draws from df_chain for each par
    a_draws = df_chain.a
    h_draws = df_chain.h
    w_draws = df_chain.w

    #scale a draws
    scaled_a_draws = a_draws .* prey_scale_adjust

    #populate matricies with posterior draws
    A_mat[i, :] = scaled_a_draws #rescale a by the prey scale adjustment factor before saving
    H_mat[i, :] = h_draws
    W_mat[i, :] = w_draws

# -------------------------------------------------------- Plotting ----------------------------------------------------------

    # Save trace plots
    plt_trace = plot(fr_chain)
    savefig(plt_trace, joinpath(trace_dir, "Inter_ID_$(i).png"))

    # Extract parameter summaries
    fitted_params = DataFrame(summarystats(fr_chain))
    fitted_a[i] = fitted_params[1,2]
    fitted_h[i] = fitted_params[2,2]
    fitted_w[i] = fitted_params[3,2]

    a_fit, h_fit, w_fit = fitted_a[i], fitted_h[i], fitted_w[i]

# -------------------------------------------------------- Empirical variance scaling plot ----------------------------------------------------

    # Plot empirical standard deviation vs average prey eaten
    plt_sigma_poff = scatter(
        empirical_sigma_plot_df.prey_offered,
        empirical_sigma_plot_df.sd_prey_eaten,
        xlabel = "Prey Offered",
        ylabel = "Empirical SD of prey eaten",
        title = "Empirical variance scaling - Prey Offered (Inter_ID $(i))",
        legend = false
    )

    plt_sigma_prdens = scatter(
        empirical_sigma_plot_df.pred_density,
        empirical_sigma_plot_df.sd_prey_eaten,
        xlabel = "Pred Density",
        ylabel = "Empirical SD of prey eaten",
        title = "Empirical variance scaling - Pred Dens (Inter_ID $(i))",
        legend = false
    )

    savefig(
        plt_sigma_poff,
        joinpath(sd_check_prey, "Empirical_SD_scaling_Prey Offered_Inter_ID_$(i).png")
    )

    savefig(
        plt_sigma_prdens,
        joinpath(sd_check_pred, "Empirical_SD_scaling_Pred_Dens_Inter_ID_$(i).png")
    )
    # -------------------------------------------------------- Surface plots ----------------------------------------------------

    # --- 3D Surface Plot (use same model used for fitting) ---
    # Establish Prey and Pred Ranges
    prey_range = range(minimum(prey_offered), stop=maximum(prey_offered), length=50)
    pred_range = range(minimum(pred_density), stop=maximum(pred_density), length=50)

    #Establish empty surface 
    Z = zeros(length(prey_range), length(pred_range))

    #Calculate posterior predcitive solution for surface based on fitting model used
    for ip in eachindex(prey_range)
        for id in eachindex(pred_range)
            preds = Float64[]

            for k in 1:length(a_draws)
                if replenished == "Y"
                    pred = predict_BDE(
                        a_draws[k],
                        h_draws[k],
                        w_draws[k],
                        prey_range[ip],
                        pred_range[id],
                        AS)
                else
                    pred = predict_BBDEL(
                        a_draws[k],
                        h_draws[k],
                        w_draws[k],
                        prey_range[ip],
                        pred_range[id],
                        AS)
                end
                push!(preds, pred)
            end

        # Posterior predictive mean
        Z[ip,id] = median(preds)

    end
end

# ----- Create Animated Surface plot -----

    # Create a surface plot 
    plt_surface = surface(
        prey_range, pred_range, Z';
        xlabel = "Prey Density",
        ylabel = "Predator Density",
        zlabel = "Prey Eaten",
        title = "Functional Response Fit (Inter_ID $(i))",
        legend = false,
        alpha = 0.7,
        zlims = (0, 1.1 * maximum(prey_eaten)),
        color = :viridis
    )

    # Overlay observed data points
    scatter3d!(plt_surface, prey_offered, pred_density, prey_eaten;
        markersize = 5, markercolor = :red, label = "Observed Data"
    )

    # Animate the rotation
    anim = @animate for angle in 0:2:360
        plot!(plt_surface, camera = (angle, 30))
    end

    # Save as GIF
    gif(anim, joinpath(animated_surface, "FR_surface_Inter_ID_$(i).gif"), fps = 15)
    
end

# ----------------------------------------------------- Save Stored Posteriors ---------------------------------------------------

# --- Save fitted parameters ---
 df_fitted = DataFrame(
    Inter_ID = 1:maximum(df_meta.Inter_ID),
    fitted_a = vec(fitted_a),
    scale_adjust_factor = prey_scale_adjust_vec,
    scaled_a = vec(fitted_a) .* prey_scale_adjust_vec,
    fitted_h = vec(fitted_h),
    fitted_w = vec(fitted_w)
)

# --- Save estimated pars in CSV ---
CSV.write("~/Inter_fitted_parameters_Revision.csv", df_fitted)

CSV.write("~/posterior_w_matrix_Revision.csv", DataFrame(W_mat, :auto))
CSV.write("~/posterior_scaled_a_matrix_Revision.csv", DataFrame(A_mat, :auto))
CSV.write("~/posterior_h_matrix_Revision.csv", DataFrame(H_mat, :auto))

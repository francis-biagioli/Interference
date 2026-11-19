function simulate_data_set_INT(loop_length,std_err,samp_size,prey_eaten,prey_offered,preds)

prey_eaten_simulated = zeros(0)
prey_offered_simulated = zeros(0)
preds_simulated = zeros(0)

for j = 1:loop_length
    std_dev_run = std_err[j].*sqrt(samp_size[j]) # calculate std dev
    # occassionally, a trial is reported as having a mean of zero and a std of 0,
    # which causes a fail on generating a distribution - assign it the meand std
    if std_dev_run == 0.0
        std_dev_run = mean(std_err).*sqrt(samp_size[j])
    end
    # make a simulation distribution from mean and std
    d = truncated(Normal(prey_eaten[j],std_dev_run),0,Inf)
    append!(prey_eaten_simulated,rand(d,Int.(samp_size[j])))
    append!(prey_offered_simulated,ones(1,Int.(samp_size[j]))*prey_offered[j])
    append!(preds_simulated,preds[j] .* ones(1,Int.(samp_size[j])))
end
prey_eaten = prey_eaten_simulated
prey_offered = prey_offered_simulated

return prey_eaten, prey_offered, preds_simulated
end

    #mu = log((prey_eaten[j]^2)/sqrt(std_dev_run^2+prey_eaten[j]^2)) # calc mu
    #sigma = sqrt(log(std_dev_run^2/(prey_eaten[j]^2)+1)) # calc sigma
    #d = LogNormal(mu,sigma)
context("Plotting methods")

test_that("plot does not crash",{
    data(NS_species_params_gears)
    data(inter)
    params <- MizerParams(NS_species_params_gears, inter)
    sim <- project(params, effort=1, t_max=4, dt = 1, t_save = 1)
    expect_that(plot(sim), not(throws_error()))
    expect_that(plotYield(sim), not(throws_error()))
    expect_that(plotYieldGear(sim), not(throws_error()))
    
    # With several resources
    params <- MizerParams(NS_species_params_gears, inter, matrix(1,12,12))
    sim <- project(params, effort=1, t_max=4, dt = 1, t_save = 1)
    expect_that(plot(sim), not(throws_error()))
})


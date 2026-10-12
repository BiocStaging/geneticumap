# create single-objective logging setup
create_logging_setup <- function(log.stats, fitness.fun.name, minimize) {
    log_stats_list <- list(fitness = list())
    for (stat in log.stats) {
        if (stat == "mean") {
            log_stats_list$fitness[[paste(fitness.fun.name, stat)]] <-
                list(fun = function(fitness) mean(as.numeric(fitness)))
        } else if (xor(stat == "worst", minimize)) {
            log_stats_list$fitness[[paste(fitness.fun.name, stat)]] <-
                list(fun = function(fitness) min(as.numeric(fitness)))
        } else {
            log_stats_list$fitness[[paste(fitness.fun.name, stat)]] <-
                list(fun = function(fitness) max(as.numeric(fitness)))
        }
    }
    return(log_stats_list)
}

create_ecr_control <- function(fitness.fun, minimize, n.objectives,
    representation, mutator, parent.selector, survival.selector) {
    control <- ecr::initECRControl(fitness.fun, minimize = minimize,
        n.objectives = n.objectives)
    control$type <- representation

    control <- ecr::registerECROperator(control, "mutate", mutator)
    control <- ecr::registerECROperator(control, "selectForMating",
        parent.selector)
    control <- ecr::registerECROperator(control, "selectForSurvival",
        survival.selector)

    return(control)
}

create_next_generation <- function(iter, control, previous_gen, lower, upper,
    lambda, p.mut, survival.strategy, n.elite) {
    parent_sel <- ecr::selectForMating(control, previous_gen$fitness,
        n.select = lambda)
    offspring <- ecr::mutate(control, previous_gen$population[parent_sel],
        p.mut = p.mut, n.iter = iter, lower = lower, upper = upper)

    # this runs in parallel batches per parallelization unit
    fitness.offspring <- ecr::evaluateFitness(control, offspring)

    for (i in seq_along(offspring)) {
        attr(offspring[[i]], "fitness") <- fitness.offspring[, i]
    }

    survival_sel <- if (survival.strategy == "plus") {
        ecr::replaceMuPlusLambda(control, previous_gen$population,
            offspring, previous_gen$fitness, fitness.offspring)
    } else {
        ecr::replaceMuCommaLambda(control, previous_gen$population, offspring,
            previous_gen$fitness, fitness.offspring, n.elite = n.elite)
    }

    # explicitly remove references to old generations
    rm(previous_gen, offspring, fitness.offspring)

    next_gen <- list(population = survival_sel$population,
        fitness = survival_sel$fitness)

    rm(survival_sel) # remove duplicate pointers

    return(next_gen)
}

evolution_loop <- function(fitness.fun, fitness.fun.name, minimize = NULL,
    n.objectives = 1L, n.dim = NULL, lower = NULL, upper = NULL, representation,
    mu, lambda, p.mut = 0.3, survival.strategy = "plus", n.elite = 0L,
    log.stats = NULL, log.pop = FALSE, monitor = NULL, initial.solutions = NULL,
    parent.selector = NULL, survival.selector = NULL, mutator = NULL,
    max.iter = 100L) {
    log_stats_list <- create_logging_setup(log.stats, fitness.fun.name,
        minimize)
    control <- create_ecr_control(fitness.fun, minimize, n.objectives,
        representation, mutator, parent.selector, survival.selector)
    log <- ecr::initLogger(control, log.stats = log_stats_list,
        log.pop = log.pop, init.size = 1000L)

    gen.fun <- ecr::genReal
    gen.pars <- list(n.dim = n.dim, lower = lower, upper = upper)

    population <- initial.solutions
    if (representation != "custom") {
        population <- do.call(ecr::initPopulation, c(list(mu = mu, gen.fun =
            gen.fun, initial.solutions = initial.solutions), gen.pars))
    }
    fitness <- ecr::evaluateFitness(control, population)

    for (i in seq_along(population)) {
        attr(population[[i]], "fitness") <- fitness[, i]
    }

    ecr::updateLogger(log, population, fitness = fitness, n.evals = mu)

    previous_gen <- list(population = population, fitness = fitness)

    for (iter in seq_len(max.iter)) {
        next_gen <- create_next_generation(iter, control, previous_gen, lower,
            upper, lambda, p.mut, survival.strategy, n.elite)

        ecr::updateLogger(log, next_gen$population, next_gen$fitness,
            n.evals = lambda)
        if (is.function(monitor)) monitor()

        previous_gen <- next_gen

        message("Finished iteration ", iter, " of ", max.iter)
    }
    return(make_ecr_result(control, log, next_gen$population, next_gen$fitness))
}

make_ecr_result <- function(control, log, population, fitness) {
    # Single-objective result object
    return(BBmisc::makeS3Obj(
        task = control$task,
        best.x = log$env$best.x,
        best.y = log$env$best.y,
        log = log,
        last.population = population,
        last.fitness = as.numeric(fitness),
        classes = c("ecr_single_objective_result", "ecr_result")
    ))
}

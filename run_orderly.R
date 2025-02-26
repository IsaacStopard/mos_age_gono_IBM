rm(list = ls())

orderly2::orderly_run("1_perennial")

orderly2::orderly_run("2_frequency")

orderly2::orderly_run("3_autocorrelation")

orderly2::orderly_run("4_process_f")

orderly2::orderly_run("5_EHT_autocorr")

orderly2::orderly_run("6_EHT_runs")

orderly2::orderly_run("7_ssm_inc_larval")

orderly2::orderly_run("8_plot")

orderly2::orderly_validate_archive(action = "orphan")
orderly2::orderly_prune_orphans()

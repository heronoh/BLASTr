condathis::create_env("bioconda::blast==1.16", env_name = "blast-env")


shell_exec <- function(cmd = "blastn", ..., env_name = "blast-env") {
  condathis::run(cmd, ..., env_name = env_name)
}

shell_exec("blastn", "-help")

temp_path <- tempfile("input_blast")

processx::run("blastn", env)

processx::process()

withr::with_envvar()

withr::local_envvar()

condathis::with_sandbox_dir()

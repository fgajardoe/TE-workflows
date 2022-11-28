ARGS=$1

CPU=6
snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.cvar.yaml Cyprinodon_variegatus.RData
snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.ochu.yaml Orestias_chungarensis.RData
snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.olau.yaml Orestias_laucaensis.RData
snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.oasc.yaml Orestias_ascotanensis.RData
snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.oglo.yaml Orestias_gloriae.RData
snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.aana.yaml Anableps_anableps.RData





ARGS=$1

CPU=6
#snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.ach.yaml Austrolebias_charrua.RData
#snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.nwhi.yaml Nematolebias_whitei.RData
#snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.cmel.yaml Cynopoecilus_melanotaenia.RData
#snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.alim.yaml Austrofundulus_limnaeus.RData
#snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.kmar.yaml Kryptolebias_marmoratus.RData
#snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.nfur.yaml Nothobranchius_furzeri.RData
#snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.olat.yaml Oryzias_latipes.RData

snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.olat.yaml panel_Aplocheilidae_genomeContrib.pdf


#snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.cvar.yaml Cyprinodon_variegatus.RData
#snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.ochu.yaml Orestias_chungarensis.RData
#snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.olau.yaml Orestias_laucaensis.RData
#snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.oasc.yaml Orestias_ascotanensis.RData
#snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.oglo.yaml Orestias_gloriae.RData
#snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.aana.yaml Anableps_anableps.RData





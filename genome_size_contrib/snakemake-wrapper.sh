ARGS=$1

CPU=6
snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.ach.yaml Austrolebias_charrua.RData
snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.nwhi.yaml Nematolebias_whitei.RData
snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.cmel.yaml Cynopoecilus_melanotaenia.RData
snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.alim.yaml Austrofundulus_limnaeus.RData
snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.nfur.yaml Nothobranchius_furzeri.RData
snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.kmar.yaml Kryptolebias_marmoratus.RData
snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.olat.yaml Oryzias_latipes.RData






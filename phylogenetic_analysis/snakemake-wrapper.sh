# borra archivos relevantes de corrida anterior
#rm *-RTE-BovB.selected.centroids.* RAxML_*all-species-RTs.RAxML 

ARGS=$1
SUFFIX=$2

CPU=6
#snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.ach.yaml ACH-RTE-BovB${SUFFIX}
snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.nwhi.yaml NWHI-RTE-BovB${SUFFIX}
snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.cmel.yaml CMEL-RTE-BovB${SUFFIX}
snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.alim.yaml ALIM-RTE-BovB${SUFFIX}
snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.nfur.yaml NFUR-RTE-BovB${SUFFIX}
#snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.kmar.yaml KMAR-RTE-BovB${SUFFIX}
snakemake --rerun-incomplete ${ARGS} -j ${CPU} -k --keep-incomplete --configfile=config.olat.yaml OLAT-RTE-BovB${SUFFIX}


# tree
#snakemake --rerun-incomplete -pj${CPU} -k --keep-incomplete --configfile=config.olat.yaml RAxML_bestTree.all-species-RTs.RAxML
#snakemake --rerun-incomplete -pj${CPU} -k --keep-incomplete --configfile=config.olat.yaml RAxML_bestTree.all-species-RTs.CDSs.RAxML

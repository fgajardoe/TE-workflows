SP=$1
ln -s ../genome_size_contrib/${SP}.gene-annotation.gff
ln -s ../genome_size_contrib/genomes/${SP}.genome.fasta
ln -s ../genome_size_contrib/genomes/${SP}.rmgr.RData

cat config.Austrolebias_charrua.yaml |perl -pe's/Austrolebias_charrua/'${SP}'/g; s/custom/ncbi/g' > config.${SP}.yaml

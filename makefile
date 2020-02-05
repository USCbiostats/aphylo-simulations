slurmr := ~/./slurmr



ALL:
	sbatch -W 01-gold-standard.sh && \
	sbatch -W 01-gold-standard-wrong-prior.sh && \
	sbatch -W 02-missinglabel.sh && \
	sbatch -W 02-missinglabel-wrong-prior.sh && \
	sbatch -W 03-misslabel.sh && \
	sbatch -W 03-misslabel-wrong-prior.sh &

bias:
	R CMD BATCH simulations/01-gold-standard/bias.r simulations/01-gold-standard/bias.Rout && \
	R CMD BATCH simulations/02-missinglabel/bias.r simulations/02-missinglabel/bias.Rout && \
	R CMD BATCH simulations/03-misslabel/bias.r simulations/03-misslabel/bias.Rout &

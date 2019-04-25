ALL:
	sbatch -W run-01-gold-standard.sh && \
	sbatch -W run-01-gold-standard-wrong-prior.sh && \
	sbatch -W run-02-missinglabel.sh && \
	sbatch -W run-02-missinglabel-wrong-prior.sh && \
	sbatch -W run-03-misslabel.sh && \
	sbatch -W run-03-misslabel-wrong-prior.sh &

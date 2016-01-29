all: shmm.so

shmm.so: shmm.cpp
	Rscript --vanilla -e 'TMB::compile("shmm.cpp")'

clean:
	rm -rf *o

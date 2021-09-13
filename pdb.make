SHELL=/usr/bin/bash

clean:
	rm -rf PDB ELS IDS BALLS BIDS

toy:
	find pdb/ -type f | head -1000 | xargs zcat | pdb load PDB IDS ELS
#	10,000 files | 18,075 models | 66,344,045 atoms | 194 types

PDB:
	find pdb/ -type f | xargs zcat | pdb load PDB IDS ELS
# 	176,191 files | 396,272 models | 1,393,323,055 atoms | 221 types

test:
	pdb test PDB IDS ELS

BALLS:
	pdb ball,r=10,CA BALLS BIDS PDB IDS ELS

PEELS:
	pdb peel PEELS BALLS

HIST:
	pdb hist,dx=0.1 HIST BINS PEELS ELS

DOCS:
	mtx DOCS = L2,idf,FS:df=2:10000000 HIST
	mtx -m 30000 transpose DOCS

LSH:
	mtx DLSH = FS:df=2:10000 HIST
	mtx LSH = simhash,lsh:L=1000,k=20 DLSH
	mtx transpose LSH

LSH.quant:
	mtx C = colsum LSH
	mtx quantiles C

transpose:
	mtx transpose DOCS
	mtx transpose LSH

LSIM:
	mtx LSIM = LSH x LSH.T top=1000

FNR: # LSH miss rate as a function of cosine similarity
	gnuplot -p -e 'FN(x,l,k) = (1-(1 - (acos(x)*2/pi))**k)**l; set grid; '\
	'p [.001:1][0:1] FN(x,100,5), FN(x,100,10), FN(x,100,20), FN(x,100,30)'

align-1:
	pdb align,near=.1,far=.5 BALLS BIDS ELS 101M:794 102M:793 | tr ' ' '\t' > z
	~/project/pdb/align.py eval < z

align-2:
	pdb dump:101M:794,hdr BALLS BIDS ELS > a1
	pdb dump:102M:793,hdr BALLS BIDS ELS > a2
	~/project/pdb/align.py align a1 a2

topK:
	~/project/pdb/align.py topK a2 > topK
	@cat topK | xcut 2 | avg.awk -vtag=1 -vH=1
	@cat topK | xcut 3 | avg.awk -vtag=2
	@cat topK | xcut 4 | avg.awk  -vtag=5
	@cat topK | xcut 5 | avg.awk -vtag=10
#     N:   mean ± std  [   min    10%    25% |  medi|    75%    90%    max] top
#     142: 1.43 ± 0.27 [  1.23   1.23   1.33 |  1.39|   1.51   1.53   3.67]   1
#     142: 1.81 ± 0.57 [  1.25   1.34   1.40 |  1.53|   2.25   2.49   4.42]   2
#     142: 2.95 ± 0.60 [  2.32   2.43   2.46 |  2.59|   3.29   3.84   5.02]   5
#     142: 3.99 ± 0.63 [  2.96   3.21   3.54 |  3.88|   4.35   4.77   6.10]  10

eval1k:
	for id in `cat REL1k.ids.2`; do \
	date '+%F,%T ------------------------------ '$$id ; \
	echo $$id | mtx QRYS2 = subset DOCS BIDS ; \
	mtx SIM2 = QRYS2 x DOCS.T top=10000 ; \
	pdb eval,top=1000,id=$$id SIM2 BALLS BIDS ELS > eval1k/$$id.atoms; \
	~/project/pdb/pdb.py eval eval1k/$$id.atoms > eval1k/$$id.eval ; \
	done

eval:
	for id in `cat REL1k.ids`; do \
	date '+%F,%T '$$id ; \
	pdb eval,top=10000,id=$$id SIM BALLS BIDS ELS > eval/$$id.atoms; \
	~/project/pdb/pdb.py eval eval/$$id.atoms > eval/$$id.eval ; \
	done

eval100k:
	for id in `cat REL100k.ids`; do \
	date '+%F,%T '$$id ; \
	pdb eval,top=100000,id=$$id SIM100k BALLS BIDS ELS > eval100k/$$id.atoms; \
	~/project/pdb/pdb.py eval eval100k/$$id.atoms > eval100k/$$id.eval ; \
	done


# relevant == distinct from query and > 50% atoms match within 0.1A
REL.rcv:
	cat eval/*.eval | gawk -vOFS='\t' '$$1 != $$2 && $$4 < 0.5 {print $$1,$$2,1}' > $@
	xcut 1 < REL.rcv | sort | uniq > REL.ids

REL1k.rcv:
	cat eval1k/*.eval | gawk -vOFS='\t' '$$1 != $$2 && $$4 < 0.5 {print $$1,$$2,1}' > $@

eplots:
	for f in eval/*.eval ; do \
	nl -v0 $$f | gawk '$$1 {print $$4, $$7}' | scatter -x11 -logx -logy ; \
	done
#

REL:
	cat REL.rcv | mtx load:rcv REL BIDS BIDS p=wrr 

QRYS:
	(cat REL.ids) | mtx QRYS = subset DOCS BIDS

SIM:
	mtx SIM = QRYS x DOCS.T top=10000

SIM100k:
	cat REL100k.ids | mtx QRYS100k = subset DOCS BIDS
	mtx SIM100k = QRYS100k x DOCS.T top=100000

log.base:
	~/project/pdb/tune-base.csh | & tee log.base

plot-base:
	grep MAP: log.base | sed 's/\.fs/ /' | sed 's/\.dx/ /' | sed 's/runs.//' | sed 's/\.evl//' \
	| gawk '{print $$1, $$7, $$8, $$9, $$2}' | matplot.py
	grep TIME: log.base | sed 's/\.fs/ /' | sed 's/\.dx/ /' | sed 's/runs.//' | sed 's/\.evl//' \
	| gawk '{print $$1, $$3, $$4, $$5, $$2}' | matplot.py

rankdistr:
	for f in runs/*.evl; do \
	gawk '{$$1=$$NF="";} 1' $$f | tr ' ' '\n' | grep . | avg.awk -vI=1 -vtag=$$f; \
	done

REL10k.evl:
	cat eval10k/*.eval \
	| gawk '$$1 != $$2 && $$5 < 0.10 {print $$1,$$2,1}' \
	| mtx load:rcv REL10k BIDS BIDS p=wrr
	mtx print:evl SIM REL10k > REL10k.evl
	cat REL10k.evl | gawk '{$$1 = $$NF = ""} 1' | tr ' ' '\n' | grep . | avg.awk -vI=1 -vH=1

REL100k.evl:
	cat eval100k/*.eval \
	| gawk '$$1 != $$2 && $$4 < 0.5 {print $$1"\t"$$2"\t1"}' \
	| mtx load:rcv REL100k BIDS BIDS p=wrr
	mtx print:evl SIM100k REL100k > REL100k.evl
	cat REL100k.evl | gawk '{$$1 = $$NF = ""} 1' | tr ' ' '\n' | grep . | avg.awk -vI=1 -vH=1

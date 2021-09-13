#!/usr/bin/tcsh -f

foreach dx (0.03 0.10 0.30 1.00)
pdb hist,dx=$dx HIST BINS PEELS ELS

foreach fs (10000 100000 1000000)

mtx TF = FS:df=2:$fs HIST

# cosine
set log = cos-tf.fs$fs.dx$dx
echo $log
mtx DOCS = weigh:L2 TF
mtx QRYS = subset DOCS BIDS < REL.ids
mtx -m 30000 transpose DOCS
/usr/bin/time -f "TIME: %e $log" mtx RET = QRYS x DOCS.T top=10000
mtx print:evl RET REL noself > runs/$log.evl
trec_eval.awk runs/$log.evl

# cosine + idf
set log = cos-idf.fs$fs.dx$dx
echo $log
mtx DOCS = weigh:idf,L2 TF
mtx QRYS = subset DOCS BIDS < REL.ids
mtx -m 30000 transpose DOCS
/usr/bin/time -f "TIME: %e $log" mtx RET = QRYS x DOCS.T top=10000
mtx print:evl RET REL noself > runs/$log.evl
trec_eval.awk runs/$log.evl

# InQuery
set log = inq.fs$fs.dx$dx
echo $log
mtx DOCS = weigh:inq TF
mtx QRYS = subset TF BIDS < REL.ids
mtx -m 30000 transpose DOCS
/usr/bin/time -f "TIME: %e $log" mtx RET = QRYS x DOCS.T top=10000
mtx print:evl RET REL noself > runs/$log.evl
trec_eval.awk runs/$log.evl

# LM (Jelinek-Mercer)
foreach b (0.01 0.03 0.10 0.30 0.90)
set log = lm$b.fs$fs.dx$dx
echo $log
mtx DOCS = weigh:lm:b=$b TF
mtx QRYS = subset TF BIDS < REL.ids
mtx -m 30000 transpose DOCS
/usr/bin/time -f "TIME: %e $log" mtx RET = QRYS x DOCS.T top=10000
mtx print:evl RET REL noself > runs/$log.evl
trec_eval.awk runs/$log.evl
end

end # fs
end # dx

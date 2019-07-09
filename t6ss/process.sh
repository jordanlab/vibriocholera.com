#!/bin/bash
set -eo pipefail
export LD_LIBRARY_PATH=/home/achande3/lib:/home/achande3/lib64:/lib:/lib6
export PATH=/home/blast/anaconda3/bin:/home/blast/anaconda3/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games
base="/home/blast/prediction_server/server/submissions/"
if [ ! -f $base/pid.lock ]; then
        echo $$ > $base/pid.lock
        i=`ls -trh  $base/ | grep ".toProcess$" -m1`
        output=`echo "${i}" | sed 's/.toProcess$//'`
        echo $i
	fasta=`ls $base/"$i" | grep ".toProcess$"`
        fname=`echo "${fasta}" | sed 's/.toProcess$//'`
        ln -sf ${base}/$i/${fasta} ${base}/$i/${fname}
        # mv "${base}/$i/$fasta" "${base}/$i/$fname"
	echo ${fasta}
        bbrename.sh in=${base}/$i/${fname} out="${base}/$i/${fasta}.fasta" prefix=sequence ignorejunk=t
        prodigal -f gff -i "${base}/$i/${fasta}.fasta" -o "${base}/$i/${fname}.gff" -a "${base}/$i/${fname}.faa" -d "${base}/$i/${fname}.fna"
        for hmm in {PAAR,aux3,lipase,vgrg,hcp,unknown,hydrolase,lysm,ntpase,transferase}; do
              hmmscan --cpu 4 -o /dev/null --tblout "${base}/$i/$hmm.output" /home/blast/prediction_server/server/hmm_profiles/${hmm}.hmm "${base}/$i/${fname}.faa"
        done
        cat  ${base}/${i}/{PAAR,aux3,lipase,vgrg,hcp,unknown,hydrolase,lysm,ntpase,transferase}.output > "${base}/$i/${fname}.output"
        echo /home/blast/prediction_server/server/parse_predictions.py -i "${base}/$i/${fname}.output" -a "${base}/$i/${fname}.faa" -n "${base}/$i/${fname}.fna" -g "${base}/$i/${fname}.gff" -f "${base}/$i/${fasta}.fasta"  -o  "${base}/$i/"
        /home/blast/prediction_server/server/parse_predictions.py -i "${base}/$i/${fname}.output" -a "${base}/$i/${fname}.faa" -n "${base}/$i/${fname}.fna" -g "${base}/$i/${fname}.gff" -f "${base}/$i/${fasta}.fasta"  -o  "${base}/$i/" | tee -a ${base}/${i}/log.txt
        cat /home/blast/prediction_server/server/done.html ${base}/${i}/*.html > /home/blast/prediction_server/server/www/results/${output}/index.html
        cp ${base}/${i}/log.txt  /home/blast/prediction_server/server/www/results/${output}/log.txt
        cp "${base}/$i/${fasta}.fasta" /home/blast/prediction_server/server/www/results/${output}/genome.fasta
        cp "${base}/$i/nucleotides.fna" /home/blast/prediction_server/server/www/results/${output}/nucleotides.fna
        cp "${base}/$i/proteins.faa" /home/blast/prediction_server/server/www/results/${output}/proteins.faa
        echo ${i}/index.html
        rm $base/pid.lock
        rm -rf "${base}/$i/"
        # /home/blast/prediction_server/server/process.sh
elif [ -f $base/pid.lock ]; then
        pid=`cat $base/pid.lock`
        if [ ! -e /proc/$pid ]; then
                rm $base/pid.lock
                /home/blast/prediction_server/server/process.sh
        fi
fi


hmmindir=$1
hmmoutdir=$2
ko=$3
name=$4
TC=$5
desc=$6

#echo $hmmindir $hmmoutdir
echo $ko $name

# name="aaeB"
# ko="K03468"
# hmmindir="./inst/extdata/hmm"
# hmmoutdir="./inst/extdata/hmm.wtc"
# TC=278.4
# desc="\"p-hydroxybenzoic acid efflux pump subunit AaeB//dbxref_kegg=K03468//dbxref_ec=NA//dbxref_TC=2.A.85.1.2\""

/usr/local/bin/gsed 's/'${ko}'/'${name}'/' ${hmmindir}'/'${ko}".hmm" > ${hmmoutdir}"/"${name}"_temp.hmm"
accession="ACC   microtraithmm:"${name}
var="TC    "${TC}" "${TC}";"
description="DESC  "${desc}
/usr/local/bin/gsed --in-place "3i$accession" ${hmmoutdir}"/"${name}"_temp.hmm"
/usr/local/bin/gsed --in-place "4i$description" ${hmmoutdir}"/"${name}"_temp.hmm"
/usr/local/bin/gsed --in-place "16i$var" ${hmmoutdir}"/"${name}"_temp.hmm"

cp ${hmmoutdir}"/"${name}"_temp.hmm" ${hmmoutdir}"/"${name}".hmm"
rm ${hmmoutdir}"/"${name}"_temp.hmm"

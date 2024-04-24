## find the different samples and the different directories 
## and concatenate the samples

## find all the sample names
cd /user/antwerpen/205/vsc20587/aitg_data/jvdabbeele/Jakke_RNAseq_immunology20240419/
find ./ -mindepth 3 -maxdepth 3  -type d |  cut -f 4 -d "/"

## get all the directories
cd /user/antwerpen/205/vsc20587/aitg_data/jvdabbeele/Jakke_RNAseq_immunology20240419/
find ${PWD}/ -mindepth 3 -maxdepth 3  -type d

## loop over the different directories
while IFS= read -r -d '' dir; do
    # Append each directory to the array
    echo ${dir}
    sample=$(basename ${dir})
    echo ${sample}

done < <(find "${PWD}/" -mindepth 3 -maxdepth 3 -type d -print0)




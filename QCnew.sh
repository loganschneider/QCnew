#!/bin/sh

# Run this QC script from within the folder of PLINK BINARY files that you want to QC
# This runs without needing a jobsub script (it's so short that it's not worth it)

module load plink
module load r

echo "Enter name of study population (e.g. WSC, MrOS, APOE), followed by [ENTER]: "
read study

echo "Enter number of genotype/chip/array pseudocohorts, followed by [ENTER]: "
read cohortnum

# gather genotype/chip/array pseudocohort(s) names into an array called "list"
for i in {1..$cohortnum}
do
	echo "Enter name(s) of genotype/chip/array pseudocohort(s), separated by spaces, followed by [ENTER]: "
	read cohortnames
	list=($cohortnames)
done

# determine phenotype name (as designated in the ${k}_pheno.txt file)
#	can generate the pheno file, using the fam2pheno.R script: https://www.dropbox.com/s/g8e5zzwkvdc8ny0/fam2pheno.R?dl=0
echo "Enter BINARY phenotype name as it appears in the {cohortname}_pheno.txt file (or leave blank), followed by [ENTER]: "
read pheno

workDir=$PWD

# loop over cohort names in "list" to submit jobs
for k in "${list[@]}"
do

#generate the long ${k}_QC.txt file, and can cat additional info as it goes
cat <<EOF >${k}_QC.txt
This is the QC log file for the ${k} cohort and ${pheno} phenotype

EOF

	#get file names based on the "cohort"/chip name, put into bfile variable
	#use that to initiate plink's first cleaning, from there the rest refer to the outputs
	#NOTE: the only files in the folder should be the origin PLINK BINARY files that you wish to QC and associated pheno files
	bfile=$(ls *.{fam,bed,bim,ped,map} -1 | sed 's/\.[a-z]*//g' | grep ${k} | head -1)
	
	#Had to add this because I believe the arbitrary or unrecorded phenotype in column 6 was making analyses confusing
	#So, I just put in the binary phenotype reported in the pheno file (assuming one is provided)
	if [ -z $pheno ]
	then
		echo "binary phenotype name not specified"
		
cat <<EOF >phenoinfam.R
famname <- list.files(pattern="^${k}.*fam$")
fam <- read.table(famname)
write.table(fam,file=paste("${k}","famWunknownphenos",sep="."),row.names=F,col.names=F,quote=F,sep=" ")
fam[,"V6"] <- -9
nonames <- unname(fam)
write.table(fam,file=famname,row.names=F,col.names=F,quote=F,sep=" ")
q()
EOF
		
		R CMD BATCH phenoinfam.R
		
	else
		echo "binary phenotype, ${pheno}, will be converted from 0/1 to 1/2 for controls/cases and then added to fam file"
		
cat <<EOF >phenoinfam.R
phenoname <- list.files(pattern="${k}_pheno.txt")
pheno <- read.table(phenoname,header=T)
write.table(pheno,file="originalphenofile.txt",row.names=F,quote=F,sep=" ")
if(2 %in% pheno[,"${pheno}"]) {
	print("case/ctrl already 2/1")
	pheno <- transform(pheno, ${pheno} = ifelse(${pheno} == 0 , -9, ${pheno}))
	pheno <- transform(pheno, ${pheno} = ifelse(${pheno} == "NA" , -9, ${pheno}))
	} else {
	pheno <- transform(pheno, ${pheno} = ifelse(${pheno} == 1 , 2, ${pheno}))
	pheno <- transform(pheno, ${pheno} = ifelse(${pheno} == 0 , 1, ${pheno}))
	pheno <- transform(pheno, ${pheno} = ifelse(${pheno} == "NA" , -9, ${pheno}))
}
write.table(pheno,file=phenoname,row.names=F,quote=F,sep=" ")
famname <- list.files(pattern="^${k}.*fam$")
fam <- read.table(famname)
write.table(fam,file=paste("${k}","famWunknownphenos",sep="."),row.names=F,col.names=F,quote=F,sep=" ")
famNOpheno <- fam[,1:5]
FIDpheno <- pheno[,c("FID","${pheno}")]
famWpheno <- merge.data.frame(famNOpheno,FIDpheno,by.x="V1",by.y="FID",sort=F)
nonames <- unname(famWpheno)
write.table(nonames,file=famname,row.names=F,col.names=F,quote=F,sep=" ")
q()
EOF
		
		R CMD BATCH phenoinfam.R
		
	fi
	
	#Need to remove SNPs with # in the name (e.g. HPA#_SNP7) as there are errors with R imports
	grep -E '#' $bfile.bim | awk '{print $2;}' > HashSNP4removal.txt #Note the limited number of SNPs with such names
	plink --bfile ${bfile} --exclude HashSNP4removal.txt --mind 0.10 --allow-no-sex --recode --out ${k}_clean_mind
	
	#grep lines of interest to be able to populate PPT cleaning flow diagram after each operation
	printf '\nBox 1: individuals removed for >10 percent missing genotype calls\n' >> ${k}_QC.txt
	grep -E 'people removed|person removed' ${k}_clean_mind.log >> ${k}_QC.txt
	sed -i -e '$a\' ${k}_QC.txt
	
	#splitting data by MAF 5 percent
	plink --file ${k}_clean_mind --maf 0.05 --allow-no-sex --recode --out ${k}_MAF_greater_5
	plink --file ${k}_clean_mind --exclude ${k}_MAF_greater_5.map --allow-no-sex --recode --out ${k}_MAF_less_5
	
	#then clean out SNPs with MAF>5 percent if missing in >5 percent of the sample and with MAF<5 percent if missing in >1 percent of samples
	plink --file ${k}_MAF_greater_5 --geno 0.05 --allow-no-sex --recode --out ${k}_MAF_greater_5_clean
	plink --file ${k}_MAF_less_5 --geno 0.01 --allow-no-sex --recode --out ${k}_MAF_less_5_clean
	
	#recombine greater and less
	plink --file ${k}_MAF_greater_5_clean --merge ${k}_MAF_less_5_clean.ped ${k}_MAF_less_5_clean.map --allow-no-sex --recode --out ${k}_MAF_clean
	
	#grep lines of interest to be able to populate PPT cleaning flow diagram after each operation
	printf '\nBox 2a: SNPs with MAF>5 percent removed if missing in >5 percent of samples\n' >> ${k}_QC.txt
	grep -E 'variants removed|variant removed' ${k}_MAF_greater_5_clean.log >> ${k}_QC.txt
	printf 'Box 2b: SNPs with MAF<5 percent removed if missing in >1 percent of samples\n' >> ${k}_QC.txt
	grep -E 'variants removed|variant removed' ${k}_MAF_less_5_clean.log >> ${k}_QC.txt
	grep -E 'people pass|person passes' ${k}_MAF_clean.log >> ${k}_QC.txt
	sed -i -e '$a\' ${k}_QC.txt
	
	#now for a stringent removal of samples missing >3 percent of their calls
	plink --file ${k}_MAF_clean --mind 0.03 --allow-no-sex --recode --out ${k}_clean2
	
	#grep lines of interest to be able to populate PPT cleaning flow diagram after each operation
	printf '\nBox 3: individuals now removed if missing >3 percent of their calls\n' >> ${k}_QC.txt
	grep -E 'people removed|person removed' ${k}_clean2.log >> ${k}_QC.txt
	grep -E 'people pass|person passes' ${k}_clean2.log >> ${k}_QC.txt
	sed -i -e '$a\' ${k}_QC.txt
	
	#check for coded and genotypic sex discrepancies
	plink --file ${k}_clean2 --check-sex --out ${k}_sex_checking
	
	#grep out the "PROBLEM" cases
	grep "PROBLEM" ${k}_sex_checking.sexcheck > ${k}_sexPROBLEM.txt
	
	#check for cryptic relatedness via identity by descent (IBD)
	rsync -aq /srv/gsfs0/projects/mignot/PLMGWAS/Files4GenAssocQCpipeline/high-ld.txt .
	#1) remove long LD-ranged regions (list of regions found here <http://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)> and put in file high-ld.txt)
	plink --file ${k}_clean2 --make-set high-ld.txt --write-set --out ${k}_hild
	plink --file ${k}_clean2 --exclude ${k}_hild.set --recode --out ${k}_trimmed
	#2) LD-prune to ensure IBD determined from SNPs only in linkage equilibrium
	plink --file ${k}_trimmed --indep 50 5 2 --out ${k}
	plink --file ${k}_trimmed --extract ${k}.prune.in --recode --out ${k}_trimmed_pruned
	###the ${k}_trimmed_pruned file is ideal for doing MDS analysis
	#3) check for cryptic relatedness
	plink --file ${k}_trimmed_pruned --genome --out ${k}_TPdups
	awk ' $10 > 0.4 ' ${k}_TPdups.genome > ${k}_TPdups.problempairs
	#4) make a text file of the likely duplicated individuals pi-hat~1
	#5) your choice as to whether you remove these individuals (i.e. what your pi-hat threshold is), but anyone with pi-hat close to 1 is likely a duplicate
	
	printf '\nBox 4: duplicate (pi-hat~1) & related (pi-hat>0.5) individuals to consider removing\n' >> ${k}_QC.txt
	cat ${k}_TPdups.problempairs >> ${k}_QC.txt
	sed -i -e '$a\' ${k}_QC.txt
	
	#Now check for excess hetero/homozygosity
	plink --file ${k}_clean2 --het --out ${k}
	
	#Generate an R script called het.R with the following
	#generates histogram of F scores
	#generates two txt files: xs homozygosity if F score is > mean+4*sd(F scores); xs heterozygosity if F score is < mean-4*sd(F scores)
cat <<EOF >het.R
Dataset <- read.table("${k}.het", header=T, sep="", na.strings="NA", dec=".", strip.white=T)
meanF<-mean(Dataset[,"F"])
sdF<-sd(Dataset[,"F"])
options(bitmapType='cairo')
png("${k}_hist.png",height=1000, width=1000)
hist(scale(Dataset[,"F"]), xlim=c(-4,4))
dev.off()
xshet <- meanF-4*sdF
xshom <- meanF+4*sdF
xshetdf <- Dataset[which(Dataset[,"F"]<xshet),]
xshomdf <- Dataset[which(Dataset[,"F"]>xshom),]
exclude <- rbind(xshetdf,xshomdf)
subexclude <- exclude[,c(1,2)]
write.table(xshetdf,file="${k}_xshet.txt",quote=F,row.names=F)
write.table(xshomdf,file="${k}_xshom.txt",quote=F,row.names=F)
write.table(subexclude,file="${k}_xshomxshet_exclude.txt",quote=F,row.names=F,col.names=F)
q()
EOF
	
	R CMD BATCH het.R
	#you can remove the individuals with xs heter/homozygosity using the ${k}_xshomxshet_exclude.txt file
	
	printf '\nBox 5: individuals to remove for xs hetero/homozygosity\n' >> ${k}_QC.txt
	printf 'XS heterozygosity (F score < mean-4*SD)\n' >> ${k}_QC.txt
	cat ${k}_xshet.txt >> ${k}_QC.txt
	sed -i -e '$a\' ${k}_QC.txt
	printf 'XS homozygosity (F score > mean+4*SD)\n' >> ${k}_QC.txt
	cat ${k}_xshom.txt >> ${k}_QC.txt
	sed -i -e '$a\' ${k}_QC.txt
	
	#next comes cleaning out CONTROLS' SNPs that are out of HWE - if a pheno file has been included
	#	can generate the pheno file, using the fam2pheno.R script: https://www.dropbox.com/s/g8e5zzwkvdc8ny0/fam2pheno.R?dl=0
	#	must the ${k}_pheno.txt file in the same folder as the files being QCed
	FILE=${k}_pheno.txt
	
	if [ -f $FILE ];
	then
		plink --file ${k}_clean2 --pheno ${k}_pheno.txt --pheno-name ${pheno} --hardy --out ${k}
		
		#Generate an R script called HWE.R with the following
cat <<EOF >HWE.R
hardy <- read.table("${k}.hwe", header=T)
names(hardy)
hwe_prob <- hardy[which(hardy[,"P"] < 1e-7),]
hwectrl <- hwe_prob[which(hwe_prob[,"TEST"]=="UNAFF"),"SNP"]
hwecase <- hwe_prob[which(hwe_prob[,"TEST"]=="AFF"),"SNP"]
write.table(hwectrl,file="${k}_ctrl_nonHWE.txt",quote=F,row.names=F,col.names=F)
write.table(hwecase,file="${k}_case_nonHWE.txt",quote=F,row.names=F,col.names=F)
q()
EOF
		
		R CMD BATCH HWE.R
		
		printf '\n' >> ${k}_QC.txt
		echo 'Box 6: SNPs removed due to departure from HWE (p<1e-7) in CONTROLS using phenotype '${pheno} >> ${k}_QC.txt
		cat ${k}_ctrl_nonHWE.txt >> ${k}_QC.txt
		sed -i -e '$a\' ${k}_QC.txt
		
		#From the whole set remove the SNPs, which are out of HWE in controls
		plink --file ${k}_clean2 --exclude ${k}_ctrl_nonHWE.txt --make-bed --out ${k}_ready4prephase
		
		#printf '\nFollowing initial genotype QC:\n' >> ${k}_QC.txt
		#grep -E 'people pass|person passes' ${k}_ready4prephase.log >> ${k}_QC.txt
	else
		printf '\n' >> ${k}_QC.txt
		echo ${FILE}' was not specified, or the phenotype is not dichotomous (preventing specification of controls)' >> ${k}_QC.txt
		sed -i -e '$a\' ${k}_QC.txt
		#printf '\nFollowing initial genotype QC:\n' >> ${k}_QC.txt
		#grep -E 'people pass|person passes' ${k}_clean2.log >> ${k}_QC.txt
		
		#Ensure uniformity of final file names for next script
		plink --file ${k}_clean2 --make-bed --out ${k}_ready4prephase
	fi
	
	#Doing the following because duplicates and other errors were found in checkVCF.py (part of prep for UMich imputation server)
	#but even these PLINK steps don't remove all of the duplicates
	plink --bfile ${k}_ready4prephase --list-duplicate-vars ids-only suppress-first --out ${k}_BPduplicates
	plink --bfile ${k}_ready4prephase --exclude ${k}_BPduplicates.dupvar --make-bed --out ${k}_NoPLINKdups
	#More BP duplicates need removal
	awk 'n=x[$1,$4]{print n"\n"$0;} {x[$1,$4]=$0;}' ${k}_NoPLINKdups.bim > ${k}_reference.dups #pulls duplicated positions by chromosome
	awk '{print $2}' ${k}_reference.dups > ${k}_Position.dups #gets a list of all rsIDs/SNP IDs for removal
	plink --bfile ${k}_NoPLINKdups --exclude ${k}_Position.dups --make-bed --out ${k}_NOdups #excludes duplicates
	
	#Include the number/ID of position duplicates removed
	cat ${k}_BPduplicates.dupvar ${k}_Position.dups | wc -l > xxx
	read IDdups < xxx
	echo ${IDdups}' rsIDs/SNP IDs were excluded, due to duplication at base positions:'  >> ${k}_QC.txt
	sed -i -e '$a\' ${k}_QC.txt
	cat ${k}_BPduplicates.dupvar >> ${k}_QC.txt
	cat ${k}_Position.dups >> ${k}_QC.txt
	sed -i -e '$a\' ${k}_QC.txt
	
	#Doing a strand orientation/alignment check
	#Used snpflip <https://github.com/biocore-ntnu/snpflip> via the SCG Cluster by <https://web.stanford.edu/group/scgpm/cgi-bin/informatics/wiki/index.php/Python>
	#Need to install it first: pip install snpflip
	# -Also housed /srv/gsfs0/projects/mignot/PLMGWAS/snpflip-master if you want it
	# -Navigate into the directory "snpflip-master"
	# -module load python
	# -pip install --user snpflip
	#OR can install in your home directory:
	# -Navigate to your home directory: /srv/gsfs0/home/{yourSUNetID}
	# -module load python/2.7
	# -mkdir -p ~/python 
	# -pip install --ignore-installed --install-option="--prefix=~/python/" /path/to/snpflip-master
	# -if doing this, the python script using snpflip must contain: module load python/2.7 export PYTHONPATH=~/python/lib/python2.7/site-packages/:$PYTHONPATH
	#Also need to download the appropriate build reference genome:
	# -from the snpflip GitHub respository's link: <http://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz>
	# -I have also stored them in the /srv/gsfs0/projects/mignot/PLMGWAS/snpflip-master directory
	
	module load python/2.7
	export PYTHONPATH=~/python/lib/python2.7/site-packages/:$PYTHONPATH
	/srv/gsfs0/home/logands/python/bin/snpflip -b $workDir/${k}_NOdups.bim -f /srv/gsfs0/projects/mignot/PLMGWAS/snpflip-master/human_g1k_v37.fasta -o $workDir/${k}_snpflip_output
	plink --bfile ${k}_NOdups --flip ${k}_snpflip_output.reverse --exclude ${k}_snpflip_output.ambiguous --make-bed --out ${k}_NOdups_aligned
	
	#Include the number/ID of reverse->flipped and ambiguous->removed SNPs
	cat ${k}_snpflip_output.reverse | wc -l > xxx
	read reversed < xxx
	echo ${reversed}' rsIDs/SNP IDs flipped to correct strand:'  >> ${k}_QC.txt
	sed -i -e '$a\' ${k}_QC.txt
	cat ${k}_snpflip_output.reverse >> ${k}_QC.txt
	sed -i -e '$a\' ${k}_QC.txt
	
	cat ${k}_snpflip_output.ambiguous | wc -l > xxx
	read ambiguous < xxx
	echo ${ambiguous}' rsIDs/SNP IDs with ambiguous strand alignment were removed:'  >> ${k}_QC.txt
	sed -i -e '$a\' ${k}_QC.txt
	cat ${k}_snpflip_output.ambiguous >> ${k}_QC.txt
	sed -i -e '$a\' ${k}_QC.txt
	
	printf '\nFollowing initial genotype QC:\n' >> ${k}_QC.txt
	grep -E 'people pass|person passes' ${k}_NOdups_aligned.log >> ${k}_QC.txt
	
done
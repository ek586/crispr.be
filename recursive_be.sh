#!/bin/bash

if [ $# -lt 2 ]; then
	echo -e "Usage:\trecursive-be [options] -v <vcf> -r <reference fa> -b <gene annotation bed>"
	exit
fi

if [ ! -d _ctmp ]; then
	mkdir _ctmp
fi

while [[ $# -ge 1 ]]; do
	key="$1"
	case $key in
		-r)
		FASTA="$2"
		shift
		;;
		-v)
		VCF="$2"
		shift
		;;
		-b)
		BED="$2"
		shift
		;;
		--default)
		;;
		*)
		;;
	esac
	shift
done

# Function to take the vcf file, outputs a gRNA design result file, and another vcf file recursively.
# The outputvcf file targets PAM candiate sites. $1: input vcf $2: output txt $3: ouput vcf
design_be()
{
# Generate SNP coordinates extended by 23 bases
	awk 'substr($1,1,1)!="#"{print $0 > "_ctmp/vcf.tmp"; print $1"\t"$2-24"\t"$2+23;}' $1 > _ctmp/region.bed

# Extract DNA sequence
	bedtools getfasta -tab -fi ${FASTA} -bed _ctmp/region.bed > _ctmp/seq.fasta

# Compare reference base with risk base, determine the base change needed, re-orient the template sequence
	paste <(cut -f2 _ctmp/seq.fasta) <(tail -n +2 $1) |
		awk -F"\t" 'function rc(seq){ rseq="";
			for(i=length(seq);i>0;--i) switch(substr(seq,i,1)) {
				case "A": rseq=rseq"T"; break;
				case "C": rseq=rseq"G"; break;
				case "G": rseq=rseq"C"; break;
				case "T": rseq=rseq"A"; break;
				default: rseq=rseq"N"; break; }
			return(rseq); }
		{$1=toupper($1);refBase=substr($1,24,1);$5=substr($5,1,1);$6=substr($6,1,1);
			if(refBase!=$6) { risk="alt"; baseChange=refBase">"$6; }
			else { risk="ref";
				switch($5) {
					case "A":
					case "C":
					case "G":
					case "T": baseChange=refBase">"$5; break;
					default:
						switch(refBase) {
							case "A": baseChange="A>G"; break;
							case "C": baseChange="C>T"; break;
							case "G": baseChange="G>A"; break;
							case "T": baseChange="T>C"; break;
							default: baseChange="N>N"; break; } }
			}
			switch(baseChange) {
				case "A>G": print $1"\t"$2"\t"$3"\t"$4"\tA>G\t"$6"\t"$9"\t+"; break;
				case "T>C": print rc($1)"\t"$2"\t"$3"\t"$4"\tT>C\t"$6"\t"$9"\t-"; break;
				case "G>A": print rc($1)"\t"$2"\t"$3"\t"$4"\tG>A\t"$6"\t"$9"\t-"; break;
				case "C>T": print $1"\t"$2"\t"$3"\t"$4"\tC>T\t"$6"\t"$9"\t+"; break;
				default: break;
			}
		}' > _ctmp/prep.txt

	# Examine PAM position sequences (+13-+17 for BE, +14-+17 for ABE) and generate output
	# Outputs a result file as $2, and a recursive vcf file as $3
	cat _ctmp/prep.txt |
		awk -F"\t" -v recprep=$3 'BEGIN{print "chr\tpos\tsnp\tchange\teditor\tgRNA\tPAM\teditinWindow\triskAllele\tInformation";
				print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > recprep; } {
			# Set PAM position range depending on BE or ABE
				if($5=="A>G"||$5=="T>C") {edit="ABE"; min_range=14;}
					else {edit="BE"; min_range=13;}
				max_range=17;
			# Prep variables
				seq=toupper($1);
				score[13]=score[14]=score[15]=score[16]=score[17]=0; score_max=1;
			# Scan through the PAM position range and score the sequences. base1 should also be G or A
				for(i=min_range;i<=max_range;++i) { pam=substr(seq,24+i,3); base1=substr(seq,4+i,1)
					if(base1 ~ /[GA]/) { if(pam ~ /.GG/) score[i]=3; else if(pam ~ /.G./) score[i]=2;
						else if(pam == "GAA" || pam == "GAT") score[i]=1;
					if(score[i]>score_max) score_max=score[i]; }}
			# Find the top score positions and save gRNA and PAM sequences
				grna=""; tpam=""; ewin="";
				for(i=min_range;i<=max_range;++i) { if(score[i]==score_max) {
					if(grna=="") { grna=substr(seq,4+i,20); tpam=substr(seq,24+i,3); ewin=21-i;}
					else { grna=grna","substr(seq,4+i,20); tpam=tpam","substr(seq,24+i,3); ewin=ewin","21-i}
				}}
			# If no position matches pam restriction, print the output first
				if(grna!="") { print $2"\t"$3"\t"$4"\t"$5"\t"edit"\t"grna"\t"tpam"\t"ewin"\t"$6"\t"$7; }
			## Find recursive editing target
			# If no position matches the original pam restriction, rescan the region
				if(grna=="") {
					for(i=min_range;i<=max_range;++i) { pam=substr(seq,24+i,3); base1=substr(seq,4+i,1)
			# First base needs to be G or A and PAM can be NAN that will be edited to NGN
						if(base1 ~ /[GA]/ && pam ~ /.A./) {
			# Report the new target pam as another vcf file and other associated information
			# On the plus strand
							if($8=="+") print $2"\t"$3+1+i"\t"$4"_relpos+"1+i"\tA\tG\t.\t.\t"$7 \
								";_;rec_gRNA"21-i";"substr(seq,4+i,20)";"pam";"edit";+" > recprep;
							else print $2"\t"$3-i-1"\t"$4"_relpos-"1+i"\tT\tC\t.\t.\t"$7 \
								";_;rec_gRNA"21-i";"substr(seq,4+i,20)";"pam";"edit";-" > recprep;
						}
					}
				}
		}' > $2
}

if true; then
design_be ${VCF} _ctmp/pres.txt _ctmp/rec1.vcf
design_be _ctmp/rec1.vcf _ctmp/pres2.txt _ctmp/rec2.vcf

# Extract cds from the gene annotation bed file
cat $BED | awk '{OFS="\t";split($11,a,","); split($12,b,","); A=""; B=""; if($7==$8) next; j=0; for(i=1;i<length(a);i++) if(($2+b[i]+a[i])>$7 && ($2+b[i])<$8) {j++; start=$2+b[i]-$7; size=a[i]; if(($2+b[i])<=$7) {start=0;size=size-($7-($2+b[i]));} if(($2+a[i]+b[i])>=$8) {size=size-($2+a[i]+b[i]-$8);} A=A""size",";B=B""start",";} print $1,$7,$8,$4,$5,$6,$7,$8,$9,j,A,B;}' > _ctmp/cds.bed

# Generate bed files of editing targets with associated information
cat _ctmp/pres2.txt | awk '\
$1!="chr" {
	split($3,a,"_relpos");
	print $1"\t"$2-1"\t"$2"\t"$3;
	if(a[2]!="") print $1"\t"$2-1-a[2]"\t"$2-a[2]"\t"a[1];
}' > _ctmp/pres.bed

# Intersect with gene annotation bed file
bedtools intersect -a _ctmp/pres.bed -b _ctmp/cds.bed -wo > _ctmp/pres.inter.tmp
bedtools intersect -a _ctmp/pres.bed -b _ctmp/cds.bed -split -wo > _ctmp/pres.interS.tmp

# Measure ths distance from the cdsStart within the block
cat _ctmp/pres.interS.tmp | awk '{ \
	exonCount = $14;
	split($15, exonSize, ",");
	split($16, exonPos, ",");
	if($10 == "+") {
		cdsStart = $6;
		cdsDist = 0;
		for(i=1; i<=exonCount; ++i) {
			if(cdsStart + exonPos[i] + exonSize[i] < $3) cdsDist += exonSize[i];
			else {
				cdsDist += $3 - (cdsStart + exonPos[i]);
				break;
			}
		}
		--cdsDist;
	} else {
		cdsEnd = $6;
		cdsDist = 0;
		for(i=exonCount; i>=1; --i) {
			if(cdsEnd + exonPos[i] > $3) cdsDist += exonSize[i];
			else {
				cdsDist += cdsEnd + exonPos[i] + exonSize[i] - $3;
				break;
			}
		}
	}
	print $0"\t"cdsDist;
}' - > _ctmp/pres.cdsDist.tmp


# Extract the bed file of the codon containing the target and get triplet sequence
cat _ctmp/pres.cdsDist.tmp | awk '{
	if($10 == "+") print $1"\t"$3-($18 % 3) -1"\t"$3-($18 % 3) + 2 "\t"$4"\t" $18":" ($18 % 3) + 1 "\t+";
	else print $1"\t"$3 + ($18 % 3) - 3"\t"$3 + ($18 % 3)"\t" $4"\t" $18 ":"($18 % 3 ) + 1 "\t-";
}' > _ctmp/pres.codon.bed

bedtools getfasta -fi ${FASTA} -bed _ctmp/pres.codon.bed -s -tab > \
	_ctmp/pres.codon.seq

# Codon table
echo "AAA	Lys	K	Lysine" > _ctmp/codon.table
echo "AAC	Asn	N	Asparagine" >> _ctmp/codon.table
echo "AAG	Lys	K	Lysine" >> _ctmp/codon.table
echo "AAT	Asn	N	Asparagine" >> _ctmp/codon.table
echo "ACA	Thr	T	Threonine" >> _ctmp/codon.table
echo "ACC	Thr	T	Threonine" >> _ctmp/codon.table
echo "ACG	Thr	T	Threonine" >> _ctmp/codon.table
echo "ACT	Thr	T	Threonine" >> _ctmp/codon.table
echo "AGA	Arg	R	Arginine" >> _ctmp/codon.table
echo "AGC	Ser	S	Serine" >> _ctmp/codon.table
echo "AGG	Arg	R	Arginine" >> _ctmp/codon.table
echo "AGT	Ser	S	Serine" >> _ctmp/codon.table
echo "ATA	Ile	I	Isoleucine" >> _ctmp/codon.table
echo "ATC	Ile	I	Isoleucine" >> _ctmp/codon.table
echo "ATG	Met	M	Methionine" >> _ctmp/codon.table
echo "ATT	Ile	I	Isoleucine" >> _ctmp/codon.table
echo "CAA	Gln	Q	Glutamine" >> _ctmp/codon.table
echo "CAC	His	H	Histidine" >> _ctmp/codon.table
echo "CAG	Gln	Q	Glutamine" >> _ctmp/codon.table
echo "CAT	His	H	Histidine" >> _ctmp/codon.table
echo "CCA	Pro	P	Proline" >> _ctmp/codon.table
echo "CCC	Pro	P	Proline" >> _ctmp/codon.table
echo "CCG	Pro	P	Proline" >> _ctmp/codon.table
echo "CCT	Pro	P	Proline" >> _ctmp/codon.table
echo "CGA	Arg	R	Arginine" >> _ctmp/codon.table
echo "CGC	Arg	R	Arginine" >> _ctmp/codon.table
echo "CGG	Arg	R	Arginine" >> _ctmp/codon.table
echo "CGT	Arg	R	Arginine" >> _ctmp/codon.table
echo "CTA	Leu	L	Leucine" >> _ctmp/codon.table
echo "CTC	Leu	L	Leucine" >> _ctmp/codon.table
echo "CTG	Leu	L	Leucine" >> _ctmp/codon.table
echo "CTT	Leu	L	Leucine" >> _ctmp/codon.table
echo "GAA	Glu	E	Glutamic_acid" >> _ctmp/codon.table
echo "GAC	Asp	D	Aspartic_acid" >> _ctmp/codon.table
echo "GAG	Glu	E	Glutamic_acid" >> _ctmp/codon.table
echo "GAT	Asp	D	Aspartic_acid" >> _ctmp/codon.table
echo "GCA	Ala	A	Alanine" >> _ctmp/codon.table
echo "GCC	Ala	A	Alanine" >> _ctmp/codon.table	
echo "GCG	Ala	A	Alanine" >> _ctmp/codon.table
echo "GCT	Ala	A	Alanine" >> _ctmp/codon.table	
echo "GGA	Gly	G	Glycine" >> _ctmp/codon.table
echo "GGC	Gly	G	Glycine" >> _ctmp/codon.table
echo "GGG	Gly	G	Glycine" >> _ctmp/codon.table
echo "GGT	Gly	G	Glycine" >> _ctmp/codon.table
echo "GTA	Val	V	Valine" >> _ctmp/codon.table
echo "GTC	Val	V	Valine" >> _ctmp/codon.table
echo "GTG	Val	V	Valine" >> _ctmp/codon.table
echo "GTT	Val	V	Valine" >> _ctmp/codon.table
echo "TAA	Stp	O	Stop" >> _ctmp/codon.table
echo "TAC	Tyr	Y	Tyrosine" >> _ctmp/codon.table
echo "TAG	Stp	O	Stop" >> _ctmp/codon.table
echo "TAT	Tyr	Y	Tyrosine" >> _ctmp/codon.table
echo "TCA	Ser	S	Serine" >> _ctmp/codon.table
echo "TCC	Ser	S	Serine" >> _ctmp/codon.table
echo "TCG	Ser	S	Serine" >> _ctmp/codon.table
echo "TCT	Ser	S	Serine" >> _ctmp/codon.table
echo "TGA	Stp	O	Stop" >> _ctmp/codon.table
echo "TGC	Cys	C	Cysteine" >> _ctmp/codon.table
echo "TGG	Trp	W	Tryptophan" >> _ctmp/codon.table
echo "TGT	Cys	C	Cysteine" >> _ctmp/codon.table
echo "TTA	Leu	L	Leucine" >> _ctmp/codon.table
echo "TTC	Phe	F	Phenylalanine" >> _ctmp/codon.table
echo "TTG	Leu	L	Leucine" >> _ctmp/codon.table
echo "TTT	Phe	F	Phenylalanine" >> _ctmp/codon.table

paste _ctmp/pres.codon.bed _ctmp/pres.codon.seq | awk 'function rc(seq) { rseq="";
	for(i=length(seq);i>0;--i) switch(substr(seq,i,1)) {
		case "A": rseq=rseq"T"; break;
		case "C": rseq=rseq"G"; break;
		case "G": rseq=rseq"C"; break;
		case "T": rseq=rseq"A"; break;
		default: rseq=rseq"N"; break; }
	return(rseq); }
NR==FNR { codon[$1] = $2; next }
{ 	split($5,a,":");
	switch(substr($8,a[2],1)) {
		case "A":
			new = substr($8,1,a[2]-1) "G" substr($8,a[2]+1,3-a[2]);
			break;
		case "C":
			new = substr($8,1,a[2]-1) "T" substr($8,a[2]+1,3-a[2]);
			break;
		case "G":
			new = substr($8,1,a[2]-1) "A" substr($8,a[2]+1,3-a[2]);
			break;
		case "T":
			new = substr($8,1,a[2]-1) "C" substr($8,a[2]+1,3-a[2]);
			break;
	}
	print $4"\t"$8"\t"codon[$8]"\t"new"\t"codon[new]; }
' _ctmp/codon.table - > _ctmp/pres.codon.res

# Merge all results and filter for candidates
awk 'NR==FNR{ n = split($1,a,"_");
	if(n==1) { tRef[a[1]] = $3; tEdt[a[1]] = $5; }
	else { pRef[a[1]] = $3; pEdt[a[1]] = $5; }
	next; }
{
	split($3,a,"_");
	if(tRef[a[1]] != tEdt[a[1]] && pRef[a[1]] && pRef[a[1]] == pEdt[a[1]]) 
		print $0"\t"tRef[a[1]]"\t"tEdt[a[1]]"\t"pRef[a[1]];
}' _ctmp/pres.codon.res _ctmp/pres2.txt | \
sort -k3,3 -u | sort -k1,1 -k2,2n > _ctmp/recRes.txt
fi

# Finalize result and ouput
awk -F"\t" 'BEGIN{
	printf "#chr\tpos\tSNPid\ttargetBaseChange\tgRNAstrand\tgRNA\tPAM\ttargetBaseEditor\ttargetBasePosition\toriginalAminoAcid\teditedAminoAcid\t";
	print "recTarget\trecBaseChange\trecStrand\trecgRNA\trecPAM\trecBaseEditor\trecBasePositioni\tPAMencodingAminoAcid\tDescription\t"}
{	split($3, a, "_relpos");
	n = split($10, b, ";");
	split($10, d, ";_;");
	split(b[n-4], c, "_gRNA");
	targetBase = substr(b[n-3],c[2],1);
	targetStrand = b[n];
	switch(targetBase) {
		case "A":
			if(targetStrand=="+") 
		 		targetBaseChange = "A>G";
			else targetBaseChange = "T>C";
		 break;
		case "C":
			if(targetStrand=="+") 
		 		targetBaseChange = "C>T";
			else targetBaseChange = "G>A";
		 break;
		case "G":
			if(targetStrand=="+") 
		 		targetBaseChange = "G>A";
			else targetBaseChange = "C>T";
		 break;
		case "T":
			if(targetStrand=="+") 
		 		targetBaseChange = "T>C";
			else targetBaseChange = "A>G";
		 break;
	}
	printf $1"\t"$2 - a[2]"\t"a[1]"\t"targetBaseChange"\t"b[n]"\t"b[n-3]"\t"b[n-2]"\t"b[n-1]"\t"c[2]"\t"$11"\t"$12"\t";
	print $2"\t"$4"\t"b[n]"\t"$6"\t"$7"\t"$5"\t"$8"\t"$13"\t"d[1];
}' _ctmp/recRes.txt
rm -rf _ctmp

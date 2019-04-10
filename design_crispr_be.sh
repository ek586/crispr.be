#!/bin/bash

if [ $# -lt 2 ]; then
	echo -e "Usage:\tdesign-crispr-be [options] -v <vcf> -r <reference fa>"
	echo -e "Options:\t-u\toutput U6 promoter compatible guides only"
	exit
fi

if [ ! -d _crispr_tmp ]; then
	mkdir _crispr_tmp
fi

U6=false

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
		-u)
		U6=true
		;;
		--default)
		;;
		*)
		;;
	esac
	shift
done

# Generate SNP coordinates extended by 23 bases
awk 'substr($1,1,1)!="#"{print $0 > "_crispr_tmp/_vcf.tmp"; print $1"\t"$2-24"\t"$2+23;}' ${VCF} > _crispr_tmp/_tmp.bed

# Extract DNA sequence
bedtools getfasta -tab -fi ${FASTA} -bed _crispr_tmp/_tmp.bed > _crispr_tmp/_tmp.fasta

# Compare reference base with risk base, determine the base change needed, re-orient the template sequence
paste <(cut -f2 _crispr_tmp/_tmp.fasta) <(tail -n +2 ${VCF}) |
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
			case "A>G": print $1"\t"$2"\t"$3"\t"$4"\tA>G\t"$6"\t"$9; break;
			case "T>C": print rc($1)"\t"$2"\t"$3"\t"$4"\tT>C\t"$6"\t"$9; break;
			case "G>A": print rc($1)"\t"$2"\t"$3"\t"$4"\tG>A\t"$6"\t"$9; break;
			case "C>T": print $1"\t"$2"\t"$3"\t"$4"\tC>T\t"$6"\t"$9; break;
			default: break;
		}
	}' > _crispr_tmp/_tmp.txt

# Examine PAM position sequences (+13-+17 for BE, +14-+17 for ABE) and generate output
cat _crispr_tmp/_tmp.txt |
	awk -F"\t" 'BEGIN{print "chr\tpos\tsnp\tchange\teditor\tgRNA\tPAM\triskAllele\tInformation"} {
		# Set PAM position range depending on BE or ABE
			if($5=="A>G"||$5=="T>C") {edit="ABE"; min_range=14;}
				else {edit="BE"; min_range=13;}
			max_range=17;
		# Prep variables
			seq=toupper($1);
			score[13]=score[14]=score[15]=score[16]=score[17]=0; score_max=1;
		# Scan through the PAM position range and score the sequences
			for(i=min_range;i<=max_range;++i) { pam=substr(seq,24+i,3);
				if(pam ~ /.GG/) score[i]=3; else if(pam ~ /.G./) score[i]=2;
					else if(pam == "GAA" || pam == "GAT") score[i]=1;
				if(score[i]>score_max) score_max=score[i]; }
		# Find the top score positions and save gRNA and PAM sequences
			grna=""; tpam="";
			for(i=min_range;i<=max_range;++i) { if(score[i]==score_max) {
				if(grna=="") { grna=substr(seq,4+i,20); tpam=substr(seq,24+i,3); }
				else { grna=grna","substr(seq,4+i,20); tpam=tpam","substr(seq,24+i,3); }
			} }
			if(grna=="") { grna="."; tpam="."; }
		# Print output
			print $2"\t"$3"\t"$4"\t"$5"\t"edit"\t"grna"\t"tpam"\t"$6"\t"$7;
	}' > _crispr_tmp/_res.txt

# Select U6
if $U6
then
	awk -F"\t" '{n=split($6,guide,",");split($7,pam,",");
		for(i=1;i<=n;++i) { firstbase=substr(guide[i],1,1);
			if(firstbase=="A"||firstbase=="G")
				print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"guide[i]"\t"pam[i]"\t"$8"\t"$9
			}
		}' _crispr_tmp/_res.txt 
else
	cat _cripsr_tmp/_res.txt
fi
rm -rf _crispr_tmp

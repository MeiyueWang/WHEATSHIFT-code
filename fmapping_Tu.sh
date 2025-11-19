function mapping(){
id=$1
genome=$2
gtf=$3
dir=$4

hisat2 -p 8 --dta -x ${genome} -1 ${dir}/${id}_tx_R1.fq.gz -2 ${dir}/${id}_tx_R2.fq.gz --summary-file ${id}.Hisat2.log |samtools view -Sb -q 20 - -o ${id}.Hisat2.q20.bam
samtools stats ${id}.Hisat2.q20.bam >${id}.samtools_stats.log

####samtools filter reads with mapping quality gt 20 ####
samtools sort -o ${id}.Hisat2.q20.sorted.bam ${id}.Hisat2.q20.bam
rm -rf ${id}.Hisat2.q20.bam
q20mappedreads=`samtools view ${id}.Hisat2.q20.sorted.bam |cut -f1 |sort -u |wc -l`
#scale=`echo ${q20mappedreads} |awk '{print 1000000/$1}'`
#bedtools genomecov -ibam ${id}.Hisat2.q20.sorted.bam -bg -split -scale ${scale} >${id}.rpm.bedgraph
#bedGraphToBigWig ${id}.rpm.bedgraph ${genome}*.fai ${id}.rpm.bw
#rm -rf ${id}.rpm.bedgraph

####featureCounts####
featureCounts -s 0 -p -t exon -a ${gtf} -o ${id}.Hisat2.featurecount.out ${id}.Hisat2.q20.sorted.bam

####stats####
cleanreads=`zcat ${dir}/${id}_tx_R1.fq.gz|wc -l |awk '{print $1/4}'`
echo -e "id\tcleanreads\tq20mappedreads" >stat_${id}.Hisat2.xls
echo -e "${id}\t${cleanreads}\t${q20mappedreads}" >>stat_${id}.Hisat2.xls
awk -v num_rawReads=${cleanreads} 'BEGIN{FS=OFS="\t"}{if(NR>1){for(i=2;i<=NF;i++){$i=$i/num_rawReads};;print $0}}' stat_${id}.Hisat2.xls >tmp_stat_${id}.xls
cat tmp_stat_${id}.xls >>stat_${id}.Hisat2.xls
rm -rf tmp_stat_${id}.xls
}

cleandir="/data/wangmeiyue/single_cell_screen/01mRNAseq/01G1812_TF_to_CS_prot/00CleanData/03batch_4"

for i in `cat /data/wangmeiyue/single_cell_screen/01mRNAseq/01G1812_TF_to_CS_prot/00CleanData/03batch_4/sampleinfo.txt |cut -f1`;do
{
mapping ${i} /data/wangmeiyue/genome/iwgsc_v1.0/hisat2Index/161010_Chinese_Spring /data/wangmeiyue/genome/iwgsc_v1.0/IWGSC_v1.1_HC_20170706_part.gtf ${cleandir}
}&
t=`expr $t + 1`
    if [ $t -ge 5 ]
    then
        echo "wait..."
        t=0
        wait
    fi
done
wait
echo "finshed!"

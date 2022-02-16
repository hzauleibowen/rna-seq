# rna-seq
Transcription analysis
183个杜洛克猪的rna-seq数据
## 1.数据下载
在linux里使用后台命令进行下载，虽然上传的wget文件里都有，但是下载的少了两个，利用R语言的setdiff函数比对原始数据和下载的数据就可以很快找到缺失的两个数据（match.R）。

```
nohup wget -O "NAME" -c site>NAME.out &  #download 
```
<img width="525" alt="微信图片_20220127211637" src="https://user-images.githubusercontent.com/82023298/151366899-41ffc129-f484-4192-aca4-25a3dae213df.png">



## 2.数据分型以及质控
在下载数据的时候发现，NCBI上有的是有fastq格式的数据的有的没有，所以统一下载NCBI格式的数据。
首先将数据名集中到一个文件夹里，按每十行来拆分成若干name文件，使用集群推荐的脚本进行处理。
```
ls S* E* >test #汇集S、E开头的文件到test文件夹里
split -| 10 test -d -a2 name_ #拆分

script.sh #集群脚本
#!/bin/bash
#SBATCH --job-name=RNA ##RNA
#SBATCH --partition=big ##作业申请的分区名称
#SBATCH --nodes=1 ##作业申请的节点数
#SBATCH --ntasks-per-node=1 ##作业申请的每个节点使用的核心数
#SBATCH --error=chaifen.err
#SBATCH --output=chaifen.out

echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"
cd /public/agis/liuyuwen_group/wangchao/double
conda init bash
conda activate RNA-seq
cat name_00 |while read id;do fastq-dump --split-3 $id;done # 拆分数据
date

sbatch scrpit.sh #提交任务
```
还有一步很重要是通过NCBI的数据来判断下载的数据是否是链特异性建库，因为链特异性建库可以区分转录本来自正义反义链的，如果不区分的化直接QC然后比对可能会造成数据的缺失。
接下来进行数据的质控，主要是去除adapter和N碱基多的reads去除低质量的片段，最后保留duplication进行下一步的定量。
在集群里面提交脚本

```
#!/bin/bash
#SBATCH --job-name=RNA ##RNA
#SBATCH --partition=low ##作业申请的分区名称
#SBATCH --nodes=1 ##作业申请的节点数
#SBATCH --ntasks-per-node=1 ##作业申请的每个节点使用的核心数
#SBATCH --error=sample_0.err
#SBATCH --output=sample_0.out
echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"
cd /public/agis/liuyuwen_group/wangchao/double/unspecificity
cat sample__0  |while read id;
do echo $id
arr=($id)
fq2=${arr[2]}
fq1=${arr[1]}
sample=${arr[0]}
trim_galore -q 25 --phred33 --length 35 -e 0.1 --stringency 4 --paired  --gzip -o clean_data  $fq1   $fq2  
done
date
```
后面就会生成质控之后的文件

![f1d49d2c764dfe8c462cfecf0ed9351](https://user-images.githubusercontent.com/82023298/153360706-c5c1e207-f29f-442e-a875-f183547a6de1.png)

下面详细记录一下脚本里面的语法
```
--quality：设定Phred quality score阈值，默认为20。
--phred33：：选择-phred33或者-phred64，表示测序平台使用的Phred quality score。
--adapter：输入adapter序列。也可以不输入，Trim Galore!会自动寻找可能性最高的平台对应的adapter。自动搜选的平台三个，也直接显式输入这三种平台，即--illumina、--nextera和--small_rna。
--stringency：设定可以忍受的前后adapter重叠的碱基数，默认为1（非常苛刻）。可以适度放宽，因为后一个adapter几乎不可能被测序仪读到。
--length：设定输出reads长度阈值，小于设定值会被抛弃。
--paired：对于双端测序结果，一对reads中，如果有一个被剔除，那么另一个会被同样抛弃，而不管是否达到标准。
--retain_unpaired：对于双端测序结果，一对reads中，如果一个read达到标准，但是对应的另一个要被抛弃，达到标准的read会被单独保存为一个文件。
--gzip和--dont_gzip：清洗后的数据zip打包或者不打包。
--output_dir：输入目录。需要提前建立目录，否则运行会报错。
-- trim-n : 移除read一端的reads

trim_galore -q 25 --phred33 --length 35 -e 0.1 --stringency 4 --paired  --gzip -o clean_data  $fq1   $fq2 #所以这个脚本大致是使用phred进行质量微调，丢弃长度为35bp的reads，与修剪序列所需的adapter序列重叠、重叠数为4，错误率小于0.1，--gzip：清洗后的数据zip打包。
```


## 3.数据比对
首先下载参考基因组
[野猪的基因组信息](http://ftp.ensembl.org/pub/release-105/fasta/sus_scrofa/dna/)
然后需要对于之前质控的文件，进行参考基因组的比对。
提交脚本
```
#!/bin/bash
#SBATCH --job-name=RNA ##RNA
#SBATCH --partition=low ##作业申请的分区名称
#SBATCH --nodes=1 ##作业申请的节点数
#SBATCH --ntasks-per-node=8 ##作业申请的每个节点使用的核心数
#SBATCH --error=sample_0.err
#SBATCH --output=sample_0.out

echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"
cd /public/agis/liuyuwen_group/wangchao/double/unspecificity/clean_data #这里是非特异性建库的质控后的数据的地址
reference=/public/agis/liuyuwen_group/wangchao/pig_refrence/Sus #参考基因组的地址
cat sample__0  |while read id;
do echo $id
arr=($id)
fq2=${arr[2]}
fq1=${arr[1]}
sample=${arr[0]}
hisat2 -x $reference -q -1 $fq1 -2 $fq2 --new-summary | samtools sort  -O bam  -@ 5 -o - > BAM/${sample}.bam  #**如果是链特异性需要加上参数** 
done
date
```
比对软件有tophat，hisat2，STAR，这里记录一下hisat2的参数
```
  hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r>} [-S <sam>]

  <ht2-idx>  Index filename prefix (minus trailing .X.ht2). 索引文件的前缀
  <m1>       Files with #1 mates, paired with files in <m2>.Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
             Paired-end测序的第一个文件，可以是压缩格式的（`.gz`或者`.bz2`）
  <m2>       Files with #2 mates, paired with files in <m1>.
             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).Paired-end测序的第二个文件
  <r>        Files with unpaired reads.
             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).single-end测序的文件
  <sam>      File for SAM output (default: stdout)
             输出比对的结果（`.sam`）

  <m1>, <m2>, <r> can be comma-separated lists (no whitespace) and can be  specified many times.  E.g. '-U file1.fq,file2.fq -U file3.fq'.可以是多个序列文件用逗号分隔

Options (defaults in parentheses):

 Input:
  -q                 query input files are FASTQ .fq/.fastq (default)输入检索序列为FASTQG格式
  --qseq             query input files are in Illumina's qseq format输入检索序列是Illumina qseq格式
  -f                 query input files are (multi-)FASTA .fa/.mfa输入文件是多个FASTA文件
  -r                 query input files are raw one-sequence-per-line输入检索序列是每行一个原始序列
  -c                 <m1>, <m2>, <r> are sequences themselves, not files直接输入序列，而不是文件
  -s/--skip <int>    skip the first <int> reads/pairs in the input (none)忽略输入序列的前`<int>`个
  -u/--upto <int>    stop after first <int> reads/pairs (no limit)只分析输入序列的前`<int>`个
  -5/--trim5 <int>   trim <int> bases from 5'/left end of reads (0)剪去reads长度为`<int>`的5'端序列
  -3/--trim3 <int>   trim <int> bases from 3'/right end of reads (0)
  
  #这里使用了smartool来将hisat2输出了sam文件转化为bam文件，samtools sort -@ 5 -o来缩小文件
  # **注意连特异性建库的时候要加上--rna-strandness RF这个参数**
```

## 4 定量
这里使用的featucount软件，RNA-seq的计数较为复杂，因为需要考虑到外显子剪切。一种计数方法是数一下与每一个被注释的外显子重合的read，另外一种方法是数一下与每一个基因区域重合的read。
提交脚本
```
#!/bin/bash
#SBATCH --job-name=RNA ##RNA
#SBATCH --partition=low ##作业申请的分区名称
#SBATCH --nodes=3 ##作业申请的节点数
#SBATCH --ntasks-per-node=8 ##作业申请的每个节点使用的核心数
#SBATCH --error=sample_1.err
#SBATCH --output=sample_1.out

echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"
cd /public/agis/liuyuwen_group/wangchao/double/unspecificity/clean_data/BAM
featureCounts -p -g gene_id -t exon -a Sus_scrofa.Sscrofa11.1.103.gtf -o featureCounts_results *.bam #p代表双端， -t表示feature-type选择外显子，-a表示参考gtf文件名，支持Gzipped文件格式，-g：当参考的gtf提供的时候，我们需要提供一个id identifier 来将feature水平的统计汇总为meta-feature水平的统计，默认为gene_id，注意！选择gtf中提供的id identifier

date
# 注释文件添加绝对路径， 特意连比对加 -s  
# 将用到的每个参数的含义写到下方
# 所有成功运行的 .err .out文件不要删除，后面需要整合一些结果
#-p 表示针对双端测序数据 ，-g指定统计的meta-feature一般是gene_id，-t指定统计的feature一般是外显子，-o输出为，-a 输出GTF\GFF的注释文件，输入前切记要添加绝对路径。
```
使用R软件对featureCount 的结果进行处理
```
#! /usr/bin/env Rscript

## Tpm
tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

## read table from featureCounts output
args <- commandArgs(T)                                            ##commandArg是R里面传递参数的函数如arg[1],arg[2]
tag <- tools::file_path_sans_ext(args[1])                         ##tools::file_path_sans_ext：r内置函数该文件包可获取不带扩展名的文件
ftr.cnt <- read.table(args[1], sep="\t", stringsAsFactors=FALSE,  ##读取表
  header=TRUE)
library(dplyr)
library(tidyr)

ftr.tpm <- ftr.cnt %>%                                            ##%>%是管道操作，来自dplyr包的管道函数，其作用是将前一步的结果直接传参给下一步的函数，从而省略了中间的赋值步骤
  gather(sample, cnt, 7:ncol(ftr.cnt)) %>%                             查看featureCount第七列是SRA开头的文件，按sample来合并，传递到定义的tmp函数里面，剔除cnt
  group_by(sample) %>%
  mutate(tpm=tpm(cnt, Length)) %>%
  select(-cnt) %>%
  spread(sample, tpm)


write.table(ftr.tpm, file=paste0(tag, "_TPM.txt"), sep="\t", row.names=FALSE, quote=FALSE)

```

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
详细记录一下trim_galore的用法
``
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


## 3.定量

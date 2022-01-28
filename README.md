# rna-seq
Transcription analysis
183个杜洛克猪的rna-seq数据
## 1.数据下载
在linux里使用后台命令进行下载，虽然上传的wget文件里都有，但是下载的少了两个，利用R语言的setdiff函数比对原始数据和下载的数据就可以很快找到缺失的两个数据（match.R）。

```
nohup wget -O "NAME" -c site>NAME.out &  #download 
```
<img width="525" alt="微信图片_20220127211637" src="https://user-images.githubusercontent.com/82023298/151366899-41ffc129-f484-4192-aca4-25a3dae213df.png">



## 2.数据分型
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
## 3.定量

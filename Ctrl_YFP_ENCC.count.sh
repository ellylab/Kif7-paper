export PATH=/home/lizhixin/softwares/cellranger-2.1.1:$PATH

sampleName=Ctrl_YFP_ENCC
appdir=/home/lizhixin/softwares/cellranger-2.1.1/
workdir=/home/lizhixin/project/scRNA-seq/rawData/10x/Aug_2018/analysis

$appdir/cellranger count --id=${sampleName}_report \
                	--transcriptome=/home/lizhixin/databases/cellranger_ref/refdata-cellranger-mm10-1.2.0 \
			--jobmode=local \
			--localcores=12 \
			--localmem=50 \
			--sample=${sampleName}_1,${sampleName}_2,${sampleName}_3,${sampleName}_4 \
			--fastqs=$workdir
			
# --csv=cellranger-tiny-bcl-simple-1.2.0.csv
# --csv=/home/lizhixin/softwares/cellranger-2.1.1/chromium-shared-sample-indexes-plate.csv \
# samplesheet

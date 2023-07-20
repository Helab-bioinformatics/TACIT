# modify T7 and T5 barcode in the "01_T7_split.sh" and "03_T5_split.sh"
mkdir 01_cut
cp ./splitbarcode0bp.pl ./01_cut
cp ./fqpick_clhuada.pl ./01_cut
cp ./01_T7_split.sh ./01_cut
cp ./02_T7_fqpick.sh ./01_cut
cp ./03_T5_split.sh ./01_cut
cp ./03_T5_split_cycle.sh ./01_cut
cp ./04_T5_fqpick.sh ./01_cut
cp ./04_T5_fqpick_cycle.sh ./01_cut

sh ./00_cut21_read1.sh
sh ./00_cut27_read2.sh

cd ./01_cut
sh 01_T7_split.sh
sh 02_T7_fqpick.sh
sh 03_T5_split_cycle.sh
sh 04_T5_fqpick_cycle.sh

mkdir h3k4me3 h3k27ac h3k27me3
find ./ -name "*raw_1_1.*" -exec cp {} ./h3k4me3 \;
find ./ -name "*raw_2_2.*" -exec cp {} ./h3k27ac \;
find ./ -name "*raw_3_3.*" -exec cp {} ./h3k27me3 \;




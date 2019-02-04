cd Mm_map
gunzip -dk GSM2183909_unique.map_95pc.txt.gz
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
do
awk -v CHR1=$CHR -v CHR2=chr$CHR 'BEGIN{FS=" ";OFS="\t"}{if($1==CHR1){$1=CHR2; print ($0)}}' \
GSM2183909_unique.map_95pc.txt > chr${CHR}_unique.map_95pc.txt
done
awk 'BEGIN{FS=" ";OFS="\t"}{if($1==20){$1="chrX"; print ($0)}}' \
GSM2183909_unique.map_95pc.txt > chrX_unique.map_95pc.txt
awk 'BEGIN{FS=" ";OFS="\t"}{if($1==21){$1="chrY"; print ($0)}}' \
GSM2183909_unique.map_95pc.txt > chrY_unique.map_95pc.txt

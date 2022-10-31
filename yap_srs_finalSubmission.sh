cd /work2/01211/cmonstr/reyap
conda activate SeqTools

# NOTE: do NOT use conda angsd here. Install the latest angsd version separately 
# (in this case angsd version is 0.938, htslib 1.15.1)


# -----remapping to chromosome-level assembly

# extracting original reads with scores from recalibrated file

# batch job:
>unbam
for F in `cat bams`;do
OUT=`echo $F | perl -pe 's/\/.+\///' | perl -pe 's/\\..+/\.fq/'`
echo $OUT
echo "samtools view $F | awk '{print \"@\"\$1\"\n\"\$10\"\n+\n\"\$11}'  >$STOCKYARD/reyap/reads/$OUT" >>unbam;
done
ls5_launcher_creator.py -j unbam -n unbam -t 0:10:00 -e matz@utexas.edu -w 24 -a tagmap -q normal
sbatch unbam.slurm


# mapping, converting to bams, indexing
>maps
REF=$WORK/db/amilV2_chroms.fasta
for F in `ls reads/*.fq`; do
OUT=`echo $F | perl -pe 's/.+\///' | perl -pe 's/\\..+//'`
echo $OUT
echo "bowtie2 --no-unal --local -x $REF -U $F -S ${OUT}.sam && samtools sort -O bam -o ${OUT}.bam ${OUT}.sam && samtools index ${OUT}.bam">>maps
done
ls5_launcher_creator.py -j maps -n maps -a tagmap -e matz@utexas.edu -t 2:00:00 -w 24 -q normal
# mapsjob=$(sbatch --dependency=afterok:$mergejob  maps.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')
mapsjob=$(sbatch maps.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')

grep overall maps.e*

# quality check
ls *[JA0-9].bam > bams.local
ls *nl.bam >bams
cat bams.local

export MinIndPerc=0.5
FILTERSQ='-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minInd $MI'
TODOQ="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"
echo 'export NIND=`cat bams.local | wc -l`; export MI=`echo "($NIND*$MinIndPerc+0.5)/1" | bc`' >calc1
echo "source calc1 && /home1/01211/cmonstr/angsd/angsd -b bams.local -r chr1:1-2000000 -GL 1 $FILTERSQ $TODOQ -P 12 -out dd && Rscript ~/bin/plotQC.R prefix=dd">a0
bash a0 &

cat quality.txt
ll -tr
wc -l bams.qc
# 291
wc -l bams
# 299

# ------ running ANGSD on each bam to get to calculate nsites per sample (and heterozygosity just because the code is this way)

>hets
>goodbams.allsites.het
REF=$STOCKYARD/db/amilV2_chroms.fasta
FILTERS='-uniqueOnly 1 -skipTriallelic 1 -minMapQ 20 -minQ 20 -doMajorMinor 5 -setMaxDepth 50'
I=0
for F in `cat bams.qc`; do
echo $I
echo "sleep $I && /home1/01211/cmonstr/angsd/angsd -i $F -anc $REF $FILTERS -GL 1 -doSaf 1 -doCounts 1 -dumpCounts 3 -out ${F/.bam/} && realSFS ${F/.bam/}.saf.idx >${F/.bam/}.ml && awk -v file=$F '{print file\"\t\"(\$1+\$2+\$3)\"\t\"\$2/(\$1+\$2+\$3)}' ${F/.bam/}.ml >>goodbams.allsites.het">>hets;
I=`echo "$I+0.5" |bc` # so lines are not overwritten
done
# check if it works (ctl-C the process if you see no errors immediately)
# head -1 hets | bash
# execute all lines in the file hets (code for TACC's stampede2)
ls6_launcher_creator.py -j hets -n hets -t 0:10:00 -e matz@utexas.edu -w 48 -a IBN21018 -q normal
sbatch hets.slurm

#---- removing high-heterozygosity outliers

R
  zz=read.table("goodbams.allsites.het")
  zz$hets=scale(zz$V3)
  zz=zz[order(zz$V3,decreasing=T),]
  zz
#                   V1      V2          V3         hets
# 49      GC_29J.bam 1512680 0.006182930  5.954880534
# 209    WOR_10A.bam 1166190 0.005957240  5.365803569
# 93    NMP_143J.bam 1132450 0.005920790  5.270664859
# 251    WOR_50A.bam 1326730 0.004485570  1.524574901
# 144    ST_118J.bam 1282010 0.004485540  1.524496598
# 62    NMP_129J.bam 1027217 0.004375644  1.237655346
# 119   NMP_159J.bam 1136522 0.004335160  1.131987430
# 243    WOR_49A.bam 1340470 0.004332200  1.124261488
# 280     WOR_9A.bam 1305730 0.004325130  1.105807972
# 222     WOR_2J.bam 1139130 0.004314890  1.079080389
# 50      GC_63A.bam 1376060 0.004291000  1.016724730
# 256    WOR_60J.bam 1225810 0.004288440  1.010042834
# 223    WOR_22J.bam 1129440 0.004287810  1.008398461
# 212  WOR_15A_1.bam 1142283 0.004286940  1.006127661
# 220     WOR_1J.bam 1252070 0.004281210  0.991171699
# 139    ST_109J.bam 1394279 0.004279980  0.987961257
# 39      GC_58J.bam 1149890 0.004277890  0.982506116
# 135    ST_117J.bam 1115070 0.004272420  0.968228784
# 236    WOR_41A.bam 1367450 0.004263010  0.943667597
# 271   WOR_6J_3.bam 1251670 0.004256780  0.927406578
# 199    ST_H762.bam 1242120 0.004255810  0.924874766
# 235    WOR_47A.bam 1065760 0.004237720  0.877657776
# 99    NMP_151J.bam 1088510 0.004236260  0.873847008
# 10     GC_106A.bam 1245046 0.004233930  0.867765439
# 31      GC_31J.bam 1132720 0.004232380  0.863719760
# 258     WOR_5J.bam 1255635 0.004229550  0.856333133
# 124    ST_112J.bam 1053750 0.004225980  0.847015021
# 44      GC_57A.bam 1200130 0.004225230  0.845057434
# 225    WOR_38A.bam  991090 0.004225210  0.845005232
# 267     WOR_8J.bam 1180990 0.004220880  0.833703432
# 246    WOR_58J.bam 1058420 0.004215450  0.819530504
# 282     WOR_8A.bam 1314848 0.004212890  0.812848609
# 18     GC_114A.bam 1248113 0.004210610  0.806897545
# 233    WOR_42A.bam 1143330 0.004209500  0.804000317
# 19     GC_120A.bam 1332383 0.004208450  0.801259696
# 27     GC_124A.bam 1084540 0.004198110  0.774271102
# 117   NMP_158J.bam 1212400 0.004195210  0.766701767
# 205    ST_H774.bam 1139003 0.004191520  0.757070441
# 101   NMP_150J.bam 1032783 0.004186410  0.743732751
# 244    WOR_52J.bam 1150470 0.004185070  0.740235196
# 128    ST_103J.bam 1283960 0.004184020  0.737494574
# 75    NMP_140J.bam  645890 0.004183640  0.736502731
# 219    WOR_20A.bam 1224380 0.004183030  0.734910560
# 64      GC_98J.bam 1153860 0.004181280  0.730342858
# 76  NMP_135J_2.bam 1157140 0.004174690  0.713142197
# 253    WOR_64J.bam 1203200 0.004173990  0.711315116
# 183    ST_H736.bam 1202280 0.004169630  0.699935013
# 113   NMP_156J.bam 1035049 0.004165100  0.688111189
# 149     ST_85J.bam 1241208 0.004164250  0.685892591
# 32      GC_61A.bam 1023000 0.004161560  0.678871381
# 224    WOR_20J.bam 1247820 0.004153950  0.659008402
# 70    NMP_134A.bam 1062440 0.004152470  0.655145431
# 8      GC_117A.bam 1059625 0.004146000  0.638257983
# 217    WOR_16J.bam 1201053 0.004144100  0.633298764
# 53      GC_66A.bam 1100573 0.004143450  0.631602189
# 13     GC_118A.bam 1223160 0.004143440  0.631576088
# 158     ST_92J.bam 1168990 0.004140530  0.623980652
# 137    ST_111J.bam 1051600 0.004138470  0.618603814
# 82    NMP_142J.bam  874806 0.004133710  0.606179664
# 136    ST_122J.bam 1099790 0.004128630  0.592920277
# 237     WOR_3J.bam 1298990 0.004128290  0.592032838
# 33      GC_32J.bam 1180315 0.004127790  0.590727780
# 226    WOR_25A.bam 1265810 0.004127790  0.590727780
# 157     ST_88J.bam 1227920 0.004123760  0.580209015
# 71    NMP_129A.bam 1327510 0.004117480  0.563817490
# 37      GC_50J.bam 1217940 0.004112210  0.550062181
# 218  WOR_15A_2.bam 1214250 0.004109810  0.543797904
# 190    ST_H737.bam 1368870 0.004108490  0.540352551
# 152     ST_89J.bam 1158610 0.004107350  0.537377020
# 69    NMP_131A.bam 1253120 0.004106200  0.534375387
# 250    WOR_65J.bam  934290 0.004101020  0.520854989
# 67    NMP_130A.bam 1106860 0.004100050  0.518323177
# 255    WOR_63J.bam 1240930 0.004097120  0.510675538
# 59    GC_70A_2.bam 1116330 0.004095940  0.507595602
# 254    WOR_54J.bam 1347170 0.004092950  0.499791357
# 28     GC_125A.bam 1211737 0.004092050  0.497442253
# 162     ST_99J.bam 1170002 0.004088950  0.489350895
# 228    WOR_22A.bam 1200470 0.004088930  0.489298692
# 264     WOR_6A.bam 1325230 0.004088630  0.488515658
# 281     WOR_9J.bam 1270040 0.004084300  0.477213858
# 111   NMP_154A.bam 1360792 0.004084150  0.476822340
# 206    WOR_14A.bam 1146800 0.004083910  0.476195912
# 227    WOR_36A.bam 1152062 0.004083600  0.475386777
# 234    WOR_24J.bam 1303000 0.004083040  0.473925112
# 118   NMP_157A.bam 1148910 0.004081190  0.469096398
# 163     ST_95J.bam 1210329 0.004081080  0.468809286
# 60    NMP_136J.bam  424393 0.004077140  0.458525431
# 133    ST_116J.bam 1147897 0.004076260  0.456228529
# 123    ST_121J.bam  928938 0.004074690  0.452130648
# 63    NMP_134J.bam  845833 0.004072580  0.446623304
# 167    ST_H729.bam 1215533 0.004070410  0.440959353
# 159   ST_93J_3.bam 1164000 0.004067880  0.434355761
# 88    NMP_140A.bam 1199460 0.004067570  0.433546625
# 213    WOR_18A.bam 1083747 0.004065020  0.426890831
# 195    ST_H745.bam 1353380 0.004058520  0.409925080
# 129    ST_102J.bam 1316340 0.004057240  0.406584132
# 55      GC_65A.bam 1275860 0.004055320  0.401572711
# 279     WOR_7J.bam 1266910 0.004053290  0.396274176
# 232    WOR_37A.bam 1260820 0.004053080  0.395726052
# 262    WOR_72J.bam 1130130 0.004052300  0.393690162
# 89    NMP_141A.bam 1140820 0.004050760  0.389670584
# 58    NMP_132J.bam  768269 0.004049900  0.387425885
# 54      GC_68A.bam 1191465 0.004047680  0.381631428
# 221    WOR_19A.bam 1268860 0.004046190  0.377742356
# 120    ST_114J.bam  965397 0.004040760  0.363569429
# 115   NMP_157J.bam 1049665 0.004039570  0.360463392
# 30      GC_27J.bam 1226450 0.004036880  0.353442181
# 86    NMP_138A.bam 1333860 0.004035550  0.349970727
# 174    ST_H732.bam 1218580 0.004035390  0.349553109
# 57    NMP_126A.bam  953600 0.004034200  0.346447071
# 22     GC_108A.bam 1277520 0.004032640  0.342375291
# 150     ST_81J.bam 1287120 0.004028840  0.332456852
# 187    ST_H748.bam 1089670 0.004024530  0.321207255
# 45      GC_60A.bam 1247670 0.004020420  0.310479680
# 131    ST_120J.bam 1021490 0.004017700  0.303380166
# 240    WOR_44A.bam 1274301 0.004017280  0.302283917
# 248    WOR_56J.bam 1268140 0.004008840  0.280254542
# 272    WOR_73J.bam 1180730 0.004008140  0.278427462
# 160    ST_H726.bam 1135208 0.004005730  0.272137083
# 108   NMP_146A.bam 1374411 0.004005640  0.271902173
# 197    ST_H764.bam 1023660 0.004004280  0.268352416
# 180    ST_H769.bam  630635 0.004003620  0.266629740
# 7      GC_109J.bam 1183290 0.004003130  0.265350783
# 238    WOR_45A.bam 1258730 0.004002870  0.264672153
# 204    WOR_13A.bam 1067338 0.004002280  0.263132185
# 1      GC_101J.bam  824776 0.003999570  0.256058772
# 125    ST_113J.bam 1114285 0.003999020  0.254623208
# 41      GC_30J.bam 1258180 0.003998000  0.251960891
# 165     ST_94J.bam 1195840 0.003995460  0.245331197
# 29      GC_33J.bam 1180510 0.003991930  0.236117489
# 109   NMP_155A.bam 1137140 0.003991820  0.235830377
# 72    NMP_132A.bam 1230050 0.003986450  0.221814057
# 77    NMP_138J.bam  850441 0.003981780  0.209624817
# 14     GC_110J.bam 1273790 0.003978620  0.201376852
# 247    WOR_55J.bam 1295030 0.003969450  0.177442093
# 186    ST_H771.bam  770618 0.003960880  0.155073403
# 90    NMP_143A.bam 1156902 0.003960730  0.154681886
# 9      GC_112A.bam 1173700 0.003959580  0.151680253
# 36      GC_48J.bam 1238960 0.003952370  0.132861321
# 4      GC_100A.bam 1086680 0.003952210  0.132443702
# 239    WOR_39A.bam 1393690 0.003951620  0.130903734
# 214  WOR_15A_3.bam 1233392 0.003951280  0.130016295
# 252    WOR_53J.bam 1372337 0.003950920  0.129076653
# 103   NMP_149A.bam 1146070 0.003949090  0.124300142
# 97    NMP_154J.bam  917175 0.003944170  0.111458374
# 65    GC_70A_1.bam 1272770 0.003940580  0.102088059
# 104   NMP_155J.bam 1041340 0.003939840  0.100156573
# 141    ST_124J.bam 1146126 0.003936270  0.090838461
# 84    NMP_142A.bam 1149181 0.003931440  0.078231603
# 26   GC_119A_1.bam 1197640 0.003929630  0.073507294
# 200    ST_H751.bam 1195310 0.003916700  0.039758501
# 105   NMP_148A.bam 1324840 0.003915010  0.035347406
# 275    WOR_71J.bam 1192360 0.003914300  0.033494224
# 107   NMP_150A.bam 1145665 0.003911420  0.025977091
# 196    ST_H747.bam 1202320 0.003898970 -0.006518847
# 276   WOR_6J_1.bam 1345130 0.003896220 -0.013696664
# 98    NMP_148J.bam 1077176 0.003893810 -0.019987043
# 15     GC_111A.bam 1211320 0.003890360 -0.028991941
# 23   GC_119A_2.bam 1172168 0.003889780 -0.030505808
# 20   GC_119A_3.bam 1175322 0.003882100 -0.050551495
# 38      GC_34J.bam 1225940 0.003881940 -0.050969114
# 173    ST_H730.bam 1229460 0.003878180 -0.060783148
# 138    ST_101J.bam 1403835 0.003875790 -0.067021324
# 178    ST_H731.bam 1166820 0.003875000 -0.069083315
# 142    ST_106J.bam 1212080 0.003873520 -0.072946286
# 155     ST_90J.bam 1196762 0.003873440 -0.073155095
# 268    WOR_75J.bam 1216830 0.003871540 -0.078114315
# 171     ST_98J.bam 1171285 0.003870460 -0.080933239
# 95    NMP_147J.bam 1021910 0.003867970 -0.087432427
# 208    ST_H775.bam 1260070 0.003866000 -0.092574355
# 188    ST_H773.bam  783083 0.003864940 -0.095341077
# 263    WOR_62J.bam 1289100 0.003862420 -0.101918568
# 47      GC_35J.bam 1302474 0.003862310 -0.102205681
# 231    WOR_25J.bam 1285460 0.003858360 -0.112515637
# 2      GC_102J.bam 1166150 0.003853790 -0.124443865
# 81    NMP_141J.bam  895804 0.003853780 -0.124469966
# 257    WOR_57J.bam 1289118 0.003848700 -0.137729353
# 168    ST_H766.bam  437118 0.003847610 -0.140574378
# 68    NMP_128A.bam 1222979 0.003840580 -0.158923490
# 241     WOR_4J.bam 1312060 0.003836330 -0.170016481
# 269    WOR_67J.bam 1332792 0.003829000 -0.189148628
# 140    ST_107J.bam 1221973 0.003828680 -0.189983865
# 106   NMP_149J.bam 1124870 0.003826040 -0.196874570
# 245    WOR_46A.bam 1386450 0.003823850 -0.202590723
# 270    WOR_70J.bam 1268399 0.003819080 -0.215040973
# 194    ST_H749.bam 1183990 0.003818670 -0.216111121
# 259    WOR_68J.bam 1145100 0.003815600 -0.224124175
# 185    ST_H738.bam 1228220 0.003810500 -0.237435764
# 151     ST_82J.bam 1224770 0.003793430 -0.281990436
# 40      GC_47J.bam 1307440 0.003788210 -0.295615239
# 6      GC_107A.bam 1246010 0.003786700 -0.299556513
# 210    WOR_12A.bam 1191060 0.003773980 -0.332757182
# 134   NMP_160J.bam 1313790 0.003771340 -0.339647887
# 143    ST_110J.bam 1394850 0.003765850 -0.353977421
# 3      GC_104A.bam 1150803 0.003765060 -0.356039412
# 145     ST_77J.bam 1215780 0.003761850 -0.364417883
# 265    WOR_69J.bam 1217465 0.003760670 -0.367497819
# 127    ST_105J.bam 1170048 0.003757810 -0.374962750
# 16     GC_121A.bam 1266540 0.003751770 -0.390727847
# 48      GC_36J.bam 1265450 0.003747620 -0.401559827
# 43    NMP_130J.bam  512591 0.003742580 -0.414714809
# 177     ST_87J.bam 1332460 0.003742540 -0.414819213
# 24     GC_113A.bam 1266580 0.003741500 -0.417533733
# 193    ST_H742.bam 1233380 0.003739800 -0.421970930
# 79    NMP_139J.bam  856558 0.003736300 -0.431106334
# 191    ST_H740.bam 1212540 0.003728860 -0.450525593
# 203    ST_H744.bam 1221370 0.003727940 -0.452926899
# 102   NMP_147A.bam 1169970 0.003727880 -0.453083506
# 35      GC_26J.bam 1315563 0.003716080 -0.483882869
# 66      GC_69A.bam 1244322 0.003714130 -0.488972594
# 202    ST_H752.bam 1227480 0.003704990 -0.512829050
# 169    ST_H727.bam 1215960 0.003700090 -0.525618616
# 184    ST_H767.bam  715780 0.003697330 -0.532822535
# 110   NMP_145A.bam 1365930 0.003688970 -0.554643100
# 201    ST_H743.bam 1376340 0.003685690 -0.563204279
# 74  NMP_135J_3.bam 1023670 0.003682780 -0.570799715
# 80    NMP_135A.bam 1197600 0.003680730 -0.576150452
# 94    NMP_146J.bam 1021523 0.003666990 -0.612013439
# 5      GC_103A.bam 1187767 0.003665690 -0.615406589
# 216    WOR_17A.bam 1045810 0.003658480 -0.634225522
# 116   NMP_158A.bam 1181140 0.003656450 -0.639524056
# 229    WOR_21A.bam 1309490 0.003653400 -0.647484908
# 198    ST_H746.bam 1325213 0.003652830 -0.648972674
# 176    ST_H733.bam 1195270 0.003649600 -0.657403347
# 92    NMP_152J.bam  791915 0.003645350 -0.668496338
# 182    ST_H734.bam 1188827 0.003640980 -0.679902543
# 166     ST_96J.bam 1199310 0.003637760 -0.688307115
# 154     ST_84J.bam 1265880 0.003632470 -0.702114626
# 260    WOR_59J.bam 1261760 0.003632030 -0.703263076
# 161     ST_91J.bam 1209740 0.003631760 -0.703967808
# 147     ST_80J.bam 1251490 0.003631310 -0.705142360
# 100   NMP_144A.bam 1181950 0.003629660 -0.709449050
# 170   ST_93J_1.bam 1274130 0.003628570 -0.712294076
# 126    ST_104J.bam 1160393 0.003626700 -0.717174992
# 12     GC_115A.bam 1136300 0.003623470 -0.725605665
# 121    ST_100J.bam 1205000 0.003622810 -0.727328341
# 25     GC_122A.bam 1147572 0.003614040 -0.750219054
# 175     ST_97J.bam 1258746 0.003596440 -0.796157087
# 21     GC_116J.bam 1269890 0.003594700 -0.800698688
# 42      GC_37J.bam 1244820 0.003594400 -0.801481722
# 122   NMP_162J.bam 1264620 0.003586750 -0.821449106
# 61      GC_67J.bam 1166060 0.003573050 -0.857207688
# 132   NMP_163J.bam 1293330 0.003567950 -0.870519277
# 112   NMP_151A.bam 1215450 0.003563220 -0.882865123
# 11     GC_105A.bam 1285243 0.003561250 -0.888007051
# 266    WOR_66J.bam 1283399 0.003555040 -0.904215868
# 146     ST_79J.bam 1239925 0.003548200 -0.922069058
# 261    WOR_61J.bam 1283477 0.003526060 -0.979857015
# 274    WOR_74J.bam 1197180 0.003513820 -1.011804829
# 181    ST_H735.bam 1201410 0.003513780 -1.011909233
# 192    ST_H741.bam 1208620 0.003493790 -1.064085442
# 52      GC_62A.bam 1326147 0.003476500 -1.109214339
# 230    WOR_24A.bam 1258677 0.003447170 -1.185769026
# 34      GC_28J.bam 1241530 0.003441550 -1.200437875
# 56    NMP_131J.bam  732914 0.003424770 -1.244235613
# 156     ST_83J.bam 1263993 0.003411310 -1.279367768
# 78    NMP_135J.bam 1159550 0.003399160 -1.311080671
# 249    WOR_48A.bam 1400006 0.003393530 -1.325775621
# 91    NMP_139A.bam 1215669 0.003392670 -1.328020321
# 207    WOR_11A.bam 1082840 0.003386680 -1.343654913
# 46      GC_49J.bam 1233770 0.003373800 -1.377273200
# 277     WOR_7A.bam 1327441 0.003371240 -1.383955096
# 83    NMP_144J.bam  761131 0.003370040 -1.387087234
# 73    NMP_137J.bam  781156 0.003369880 -1.387504853
# 153     ST_86J.bam 1243560 0.003358970 -1.415981213
# 278   WOR_6J_2.bam 1326210 0.003352940 -1.431720209
# 164   ST_93J_2.bam 1177650 0.003343430 -1.456542408
# 242    WOR_43A.bam 1275234 0.003339000 -1.468105219
# 85    NMP_137A.bam 1255950 0.003282120 -1.616568588
# 148     ST_78J.bam 1191190 0.003251950 -1.695315773
# 211    WOR_16A.bam  951702 0.003251740 -1.695863897
# 189    ST_H739.bam 1169500 0.003232990 -1.744803563
# 179    ST_H728.bam 1169894 0.003232660 -1.745664901
# 51      GC_56A.bam 1295510 0.003113310 -2.057182185
# 215    WOR_15A.bam 1067190 0.003074520 -2.158428565
# 87    NMP_145J.bam  837293 0.003051280 -2.219087649
# 273     WOR_6J.bam 1095540 0.003034800 -2.262102352
# 96    NMP_153J.bam  774481 0.002995230 -2.365384622
# 114   NMP_156A.bam 1119420 0.002929630 -2.536608199
# 172     ST_93J.bam  906970 0.002854390 -2.732993288
# 130   NMP_161J.bam 1089824 0.002696720 -3.144530198
# 17     GC_119A.bam  840977 0.002540240 -3.552961070

  # eyeball where the break is. In this case, there are three obvious high-het outliers more than 5 SDs higher.
  highH=zz$V1[1:3]
  bams=scan("bams.qc",what="c")
  bams=bams[!(bams %in% highH)]
  length(bams)
  write(bams,"goodbams.nohh")
  q()
  
#----------------------- quick and dirty IBS for detecting clones and wrong collections, with minInd=[75% of total]

idev
FILTERS='-minInd 216 -uniqueOnly 1 -skipTriallelic 1 -minMapQ 20 -minQ 20 -doHWE 1 -maxHetFreq 0.5 -snp_pval 1e-5 -sb_pval 1e-3 -hetbias_pval 1e-3'
TODO='-doMajorMinor 1 -dosnpstat 1 -doMaf 1 -doCounts 1 -doPost 1 -doGlf 3 -makeMatrix 1 -doIBS 1'
echo "/home1/01211/cmonstr/angsd/angsd -b goodbams.nohh -GL 1 $FILTERS $TODO -P 8 -out all.nohh">all.nohh
bash all.nohh
        # -> Total number of sites analyzed: 3312932
        # -> Number of sites retained after filtering: 41216

#----------- detecting and removing clones and wrong collections

# copy all.nohh.ibsMat and goodbams.nohh to laptop,run ibs_all_removeBads.R on it interactively, copy bams2 back to HPC

#---- which of the remaining samples contain the least sites?

R
  zz=read.table("goodbams.allsites.het")
  bams=scan("bams2",what="c")
  zz=zz[zz$V1 %in% bams,]
  zz=zz[order(zz$V2),]
  head(zz)
  write(zz$V1[1:3],"3worst")
  q()

# ------- extracting sites shared between three lowest-sites samples

REF=$STOCKYARD/db/amilV2_chroms.fasta
FILTERS='-minInd 3 -uniqueOnly 1 -skipTriallelic 1 -minMapQ 20 -minQ 20 -doHWE 1 -maxHetFreq 0.5 -sb_pval 1e-3 -hetbias_pval 1e-3'
TODO="-doSaf 1 -doCounts 1 -anc $REF -ref $REF -doMajorMinor 1 -doMaf 1 -doGeno 1 -dosnpstat 1 -doPost 1 -doGlf 2 -makeMatrix 1 -doIBS 1"
echo "/home1/01211/cmonstr/angsd/angsd -b 3worst -GL 1 -P 6 $FILTERS $TODO -out 3w">3w
bash 3w
        # -> Total number of sites analyzed: 702634
        # -> Number of sites retained after filtering: 257264

# saving and indexing sites
zcat 3w.mafs.gz | cut -f 1,2 | tail -n +2 >goodsites
/home1/01211/cmonstr/angsd/angsd sites index goodsites

#----------------------- for only good bams (no clones and wrong collections)

# glf3 for ngsRelate and IBS 
FILTERS='-minInd 195 -uniqueOnly 1 -skipTriallelic 1 -minMapQ 20 -minQ 20 -doHWE 1 -maxHetFreq 0.5 -snp_pval 1e-5 -sb_pval 1e-3 -hetbias_pval 1e-3'
TODO1='-doMajorMinor 1 -dosnpstat 1 -doMaf 1 -doCounts 1 -doPost 1 -doGlf 3 -makeMatrix 1 -doIBS 1'
echo "/home1/01211/cmonstr/angsd/angsd -b bams2 -sites goodsites -GL 1 $FILTERS $TODO1 -P 8 -out g3">g3
bash g3
        # -> Total number of sites analyzed: 3222254
        # -> Number of sites retained after filtering: 11089

# ----- relatedness with NgsRelate

zcat g3.mafs.gz | cut -f5 |sed 1d >freq
ngsRelate -g g3.glf.gz -n 261 -f freq -O g3.relatedness

# ----- pcangsd

# rerunning with -dpGlf 2 (beagle output)
FILTERS='-minInd 195 -uniqueOnly 1 -skipTriallelic 1 -minMapQ 20 -minQ 20 -doHWE 1 -maxHetFreq 0.5 -snp_pval 1e-5 -sb_pval 1e-3 -hetbias_pval 1e-3'
TODO1='-doMajorMinor 1 -dosnpstat 1 -doMaf 1 -doCounts 1 -doPost 1 -doGlf 2 '
echo "/home1/01211/cmonstr/angsd/angsd -b bams2 -sites goodsites -GL 1 $FILTERS $TODO1 -P 8 -out g2">g2
bash g2 &&\
pcangsd --beagle g2.beagle.gz --admix -o pcangsd.min --inbreedSamples -t 12 --minMaf 0.025 --selection --sites_save --maf_save &

# ---- heterozygosity on min sites only

>hets
REF=$STOCKYARD/db/amilV2_chroms.fasta
FILTERS='-uniqueOnly 1 -skipTriallelic 1 -minMapQ 20 -minQ 20 -doMajorMinor 5 -setMaxDepth 50'
ANGSDPATH='/home1/01211/cmonstr/angsd'
I=0
for F in `cat bams2`; do
echo "${ANGSDPATH}/angsd -sites goodsites -i $F -anc $REF $FILTERS -GL 1 -doSaf 1 -doCounts 1 -dumpCounts 3 -out ${F/.bam/} && realSFS ${F/.bam/}.saf.idx >${F/.bam/}.ml && awk -v file=$F '{print file\"\t\"(\$1+\$2+\$3)\"\t\"\$2/(\$1+\$2+\$3)}' ${F/.bam/}.ml >${F/.bam/}.goodsites.het">>hets;
done
# check if it works (ctl-C the process if you see no errors immediately)
# head -1 hets | bash
# execute all lines in the file hets (code for TACC's stampede2)
ls6_launcher_creator.py -j hets -n hets -t 0:10:00 -e matz@utexas.edu -w 48 -a IBN21018 -q normal
sbatch hets.slurm
cat *.goodsites.het >goodsites.het

#---- calculating thetas per chromosome for each pop_age group

REF=/work/01211/cmonstr/db/amilV2_chroms.fasta
TODO="-doSaf 1 -anc $REF -ref $REF"
ANGSDPATH='/home1/01211/cmonstr/angsd'
# echo "chr1
# chr2">c12

# requires winsfs, https://github.com/malthesr/winsfs
>sfsj
for F in `ls *.popbams`; do
echo "${ANGSDPATH}/angsd -b $F -sites goodsites -GL 1 -P 8 $FILTERS $TODO -out $F && \
${ANGSDPATH}/misc/realSFS ${F}.saf.idx > ${F}.sfs && ${ANGSDPATH}/misc/realSFS saf2theta ${F}.saf.idx -outname $F -sfs ${F}.sfs && \
${ANGSDPATH}/misc/thetaStat do_stat ${F}.thetas.idx -outnames $F && grep -v \"WinCenter\" ${F}.pestPG | awk '{ print \$5/\$14}' >${F}_piPerChrom">>sfsj;
done
ls6_launcher_creator.py -j sfsj -n sfsj -t 0:20:00 -e matz@utexas.edu -w 4 -a IBN21018
sbatch sfsj.slurm
cat sfsj

# rendomized bam lists
REF=/work/01211/cmonstr/db/amilV2_chroms.fasta
TODO="-doSaf 1 -anc $REF -ref $REF"
ANGSDPATH='/home1/01211/cmonstr/angsd'
>sfsj
for F in `ls *.popbamsR`; do
echo "${ANGSDPATH}/angsd -b $F -sites goodsites -GL 1 -P 8 $FILTERS $TODO -out $F && \
${ANGSDPATH}/misc/realSFS ${F}.saf.idx > ${F}.sfs && ${ANGSDPATH}/misc/realSFS saf2theta ${F}.saf.idx -outname $F -sfs ${F}.sfs && \
${ANGSDPATH}/misc/thetaStat do_stat ${F}.thetas.idx -outnames $F && grep -v \"WinCenter\" ${F}.pestPG | awk '{ print \$5/\$14}' >${F}_piPerChrom">>sfsj;
done
ls6_launcher_creator.py -j sfsj -n sfsj -t 0:20:00 -e matz@utexas.edu -w 4 -a IBN21018
sbatch sfsj.slurm
head -1 sfsj | bash


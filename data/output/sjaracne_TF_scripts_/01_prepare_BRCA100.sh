perl -pe 's/\r\n|\n|\r/\n/g' data/tf.txt > data/tf.txt.tmp
rm data/tf.txt
mv data/tf.txt.tmp data/tf.txt
perl -pe 's/\r\n|\n|\r/\n/g' data/BRCA100.exp > data/BRCA100.exp.tmp
rm data/BRCA100.exp
mv data/BRCA100.exp.tmp data/BRCA100.exp

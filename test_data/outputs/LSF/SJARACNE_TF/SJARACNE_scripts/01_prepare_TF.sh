perl -pe 's/\r\n|\n|\r/\n/g' ./test_data/inputs/tf.txt > ./test_data/inputs/tf.txt.tmp
rm ./test_data/inputs/tf.txt
mv ./test_data/inputs/tf.txt.tmp ./test_data/inputs/tf.txt
perl -pe 's/\r\n|\n|\r/\n/g' ./test_data/inputs/BRCA100.exp > ./test_data/inputs/BRCA100.exp.tmp
rm ./test_data/inputs/BRCA100.exp
mv ./test_data/inputs/BRCA100.exp.tmp ./test_data/inputs/BRCA100.exp
